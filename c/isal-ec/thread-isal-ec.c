#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <stdio.h>
#include <limits.h>
#include <libaio.h>
#include <stdint.h>
#include <sys/time.h>
#include <utime.h>
#include <errno.h>
#include <isa-l.h>
#include <pthread.h>
#include <sys/sysinfo.h>

#include "common.h"
#include "error.h"

#define M_K_P_MAX	255
#define K_DEFAULT	2
#define P_DEFAULT	1

typedef unsigned char u8;

typedef struct erasure_code_buf_info {
	int		m;
	int		k;
	int		p;
	int		frag_len;
	/* ec buffer */
      	unsigned char	*frag_ptrs[M_K_P_MAX];
	unsigned char	*recover_srcs[M_K_P_MAX];
	unsigned char	*recover_outp[M_K_P_MAX];
        unsigned char	frag_err_list[M_K_P_MAX];
	int		nerrs;
      	unsigned char	*recover_frag_ptrs[M_K_P_MAX];
	
	// Coefficient matrices
	unsigned char	*encode_matrix;
	unsigned char	*decode_matrix;
	unsigned char 	*invert_matrix;
	unsigned char	*temp_matrix;
	unsigned char	*g_tbls;
	unsigned char	decode_index[M_K_P_MAX];
} EC_BUF_INFO;

EC_BUF_INFO* alloc_ec_buf(int m, int k, int p, int frag_len)
{
	int	i;
	EC_BUF_INFO	*ebi;

	ebi = malloc(sizeof(EC_BUF_INFO));
	if (ebi == NULL) {
		ERR_SYS("malloc() error");
	}
	memset(ebi, 0, sizeof(EC_BUF_INFO));
	ebi->m = m;
	ebi->k = k;
	ebi->p = p;
	ebi->frag_len = frag_len;
	ebi->nerrs = 0;
	
        // Allocate coding matrices
        ebi->encode_matrix = malloc(m * k);
        ebi->decode_matrix = malloc(m * k);
        ebi->invert_matrix = malloc(m * k);
        ebi->temp_matrix = malloc(m * k);
        ebi->g_tbls = malloc(k * p * 32);

        if (ebi->encode_matrix == NULL || ebi->decode_matrix == NULL || ebi->invert_matrix == NULL
		|| ebi->temp_matrix == NULL || ebi->g_tbls == NULL) {
		ERR_SYS("malloc() error");
        }
	// alloc for frag
	for (i = 0; i < m; i++) {
		ebi->frag_ptrs[i] = malloc(ebi->frag_len);
		if (ebi->frag_ptrs[i] == NULL) {
			ERR_SYS("malloc() error");
		}
	}
	// alloc for recover_outp
	for (i = 0; i < p; i++) {
		ebi->recover_outp[i] = malloc(ebi->frag_len);
		if (ebi->recover_outp[i] == NULL) {
			ERR_SYS("malloc() error");
		}
	}
	return ebi;
}

void release_ec_buf(EC_BUF_INFO *ebi)
{
	int	i;

	if (ebi == NULL) {
		return;
	}
	if (ebi->encode_matrix) free(ebi->encode_matrix);
	if (ebi->decode_matrix) free(ebi->decode_matrix);
	if (ebi->invert_matrix) free(ebi->invert_matrix);
	if (ebi->temp_matrix) free(ebi->temp_matrix);
	if (ebi->g_tbls) free(ebi->g_tbls);

	for (i = 0; i < ebi->m; i++) {
		if (ebi->frag_ptrs[i] != NULL) {
			free(ebi->frag_ptrs[i]);
		}
	}
	for (i = 0; i < ebi->p; i++) {
		if (ebi->recover_outp[i] != NULL) {
			free(ebi->recover_outp[i]);
		}
	}
	free(ebi);

	//dbg("release all buffer");
}


static int gf_gen_decode_matrix_simple(u8 * encode_matrix,
                                       u8 * decode_matrix,
                                       u8 * invert_matrix,
                                       u8 * temp_matrix,
                                       u8 * decode_index, u8 * frag_err_list, int nerrs, int k,
                                       int m)
{
        int i, j, p, r;
        int nsrcerrs = 0;
        u8 s, *b = temp_matrix;
        u8 frag_in_err[M_K_P_MAX];

        memset(frag_in_err, 0, sizeof(frag_in_err));

        // Order the fragments in erasure for easier sorting
        for (i = 0; i < nerrs; i++) {
                if (frag_err_list[i] < k)
                        nsrcerrs++;
                frag_in_err[frag_err_list[i]] = 1;
        }

        // Construct b (matrix that encoded remaining frags) by removing erased rows
        for (i = 0, r = 0; i < k; i++, r++) {
                while (frag_in_err[r])
                        r++;
                for (j = 0; j < k; j++)
                        b[k * i + j] = encode_matrix[k * r + j];
                decode_index[i] = r;
        }

        // Invert matrix to get recovery matrix
        if (gf_invert_matrix(b, invert_matrix, k) < 0)
                return -1;

        // Get decode matrix with only wanted recovery rows
        for (i = 0; i < nerrs; i++) {
                if (frag_err_list[i] < k)       // A src err
                        for (j = 0; j < k; j++)
                                decode_matrix[k * i + j] =
                                    invert_matrix[k * frag_err_list[i] + j];
        }

        // For non-src (parity) erasures need to multiply encode matrix * invert
        for (p = 0; p < nerrs; p++) {
                if (frag_err_list[p] >= k) {    // A parity err
                        for (i = 0; i < k; i++) {
                                s = 0;
                                for (j = 0; j < k; j++)
                                        s ^= gf_mul(invert_matrix[j * k + i],
                                                    encode_matrix[k * frag_err_list[p] + j]);
                                decode_matrix[k * p + i] = s;
                        }
                }
        }
        return 0;
}

int time_since(struct timeval *tvold)
{
	struct timeval	tvnow;
	int		gap;

	gettimeofday(&tvnow, NULL);
	gap = (tvnow.tv_sec - tvold->tv_sec)*1000*1000 + (tvnow.tv_usec - tvold->tv_usec);

	return gap;	// return 'us'
}

typedef struct erasure_sharding_index_info {
        int     m;
        int     k;
        int     p;
        int     frag_len;
        int     fd;
	int	wfd[M_K_P_MAX];
	int	block_len;
	int64_t	offset;
        int     index;
        int	time;
} THREAD_BLOCK_INFO;

int		nr_cpus;
void *pthread_encode_ec_block(void *arg)
{
	THREAD_BLOCK_INFO *t_block_info = (THREAD_BLOCK_INFO *)arg;
	int	m, k, p, frag_len, index;
	int	fd, i;
	int64_t	offset;
	EC_BUF_INFO	*ebi;
	struct timeval	start;
	
	m = t_block_info->m;
	k = t_block_info->k;
	p = t_block_info->p;
	frag_len = t_block_info->frag_len;
	index = t_block_info->index;
	fd = t_block_info->fd;
	offset = t_block_info->offset;


	DBG("ptid[%ld], m[%d], k[%d], p[%d], frag_len[%d], offset[%lld], index[%d]", pthread_self(), m, k, p, frag_len, offset, index);
	ebi = alloc_ec_buf(m, k, p, frag_len);
	for (i = 0; i < k; i++) {
		if (preadn(fd, ebi->frag_ptrs[i], frag_len, offset + i * frag_len) != frag_len) {
			ERR_SYS("preadn() error)");
		}
	}

	gettimeofday(&start, NULL);
	gf_gen_cauchy1_matrix(ebi->encode_matrix, m, k);
	// Initialize g_tbls from encode matrix
	ec_init_tables(k, p, &(ebi->encode_matrix)[k * k], ebi->g_tbls);
	// Generate EC parity blocks from sources
	ec_encode_data(frag_len, k, p, ebi->g_tbls, ebi->frag_ptrs, &(ebi->frag_ptrs)[k]);

	t_block_info->time = time_since(&start);
	for (i = 0; i < ebi->m; i++) {
		DBG("ptid[%ld], wfd[%d]: %d, index[%d], pwriten(offset): %lld, write_len:[%ld]", pthread_self(), i, t_block_info->wfd[i], index, index*frag_len, ebi->frag_len);
		if(pwriten(t_block_info->wfd[i], ebi->frag_ptrs[i], ebi->frag_len, index * frag_len) != ebi->frag_len) {
			ERR_SYS("writen() error");
		}
	}
	release_ec_buf(ebi);
	
	return NULL;
}

int main(int argc, char **argv)
{
	int             opt;
	int		fd, tmpfd, wfd[M_K_P_MAX] = {0,};
	struct stat	st;
	int64_t		file_size, block_len, frag_len;
	int		m, k, p;
	int		i, is_decode, total_time;
	char		filename[NAME_MAX], tmpname[NAME_MAX];
	EC_BUF_INFO	*ebi;
	struct timeval	start;
	pthread_t 	*ptid;
	THREAD_BLOCK_INFO	*t_block_info;

	is_decode = 0;
	k = K_DEFAULT;
	p = P_DEFAULT;
        while ((opt = getopt(argc, argv, "dk:p:")) != -1)
        {
                switch (opt)
                {
                        case 'd':
                                is_decode = 1;
                                break;
                        case 'k':
                                k = strtoul(optarg, NULL, 10);
                                break;
                        case 'p':
                                p = strtoul(optarg, NULL, 10);
                                break;
                        default:
                		err_quit("USAGE: %s [-d] [-k k] [-p p] <origin_file | encode_file_prefix>", argv[0]);
				break;
                }
        }
	m = k + p;
	if (m >= M_K_P_MAX || k < 1 || p < 1 ) {
		err_quit("invalid parameters: (k+p)[%d] or k [%d] or p[%d] invalid", m, k, p);
	}
	if (argc - optind != 1) {
                err_quit("USAGE: %s [-d] [-k k] [-p p] <origin_file | encode_file_prefix>", argv[0]);
        }

	nr_cpus = get_nprocs();
	ptid = malloc(nr_cpus * sizeof(pthread_t));
	if (NULL == ptid) {
		ERR_SYS("malloc(pthread_t) error");
	}
	memset(ptid, 0, nr_cpus * sizeof(pthread_t));
	t_block_info = malloc(nr_cpus * sizeof(THREAD_BLOCK_INFO));
	if (NULL == t_block_info) {
		ERR_SYS("malloc(THREAD_BLOCK_INFO) error");
	}
	memset(t_block_info, 0, nr_cpus * sizeof(THREAD_BLOCK_INFO));

	file_size = 0;
	frag_len = 0;
	strncpy(filename, argv[optind], sizeof(filename));
	if (is_decode == 0) {
		if (lstat(filename, &st) < 0) {
			ERR_SYS("lstat('%s') error", filename);
		}
		if ((fd = open(filename, O_RDONLY)) < 0) {
			ERR_SYS("open('%s') error", filename);
		}
		file_size = st.st_size;
		block_len = file_size / nr_cpus;
		frag_len = block_len / k;
		//msg("file['%s'], file_size[%ld], m[%d], k[%d], p[%d], block_len[%ld], frag_len[%ld], nr_cpus[%d]", filename, file_size, m, k, p, block_len, frag_len, nr_cpus);

		for (i = 0; i < m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.%d", filename, i);
			if ((wfd[i] = open(tmpname, O_CREAT | O_RDWR, 0644)) < 0) {
				ERR_SYS("open('%s') error", tmpname);
			}
			DBG("file:[%s], wfd[%d]: %d", tmpname, i, wfd[i]);
		}

		for (i = 0; i < nr_cpus; i++) {
			t_block_info[i].m = m;
			t_block_info[i].k = k;
			t_block_info[i].p = p;
			t_block_info[i].frag_len = frag_len;
			t_block_info[i].block_len = block_len;
			t_block_info[i].fd = fd;
			memcpy(t_block_info[i].wfd, wfd, sizeof(wfd));
			t_block_info[i].offset = i*block_len;
			t_block_info[i].index = i;
			t_block_info[i].time = 0;

			pthread_create(&ptid[i], NULL, pthread_encode_ec_block, &t_block_info[i]);
		}	

		for (i = 0; i < nr_cpus; i++) {
			pthread_join(ptid[i], NULL);
		}	
		close(fd);
		for ( i = 0; i < m; i++) {
			close(wfd[i]);
		}
		total_time = 0;
		for (i = 0; i < nr_cpus; i++) {
			if (t_block_info[i].time > total_time) {
				total_time = t_block_info[i].time;
			}
		}
		msg("COST TIME: %lld (us)", total_time);
	}
	else {
		// get frag_len
		for (i = 0; i < m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.%d", filename, i);
			if (lstat(tmpname, &st) < 0) {
				DBG("lstat('%s') error, skip it!", tmpname);
				continue;
			}
			frag_len = st.st_size;
			break;
		}
		ebi = alloc_ec_buf(m, k, p, frag_len);
		dbg("m[%d], k[%d], p[%d], frag_len[%d]", ebi->m, ebi->k, ebi->p, ebi->frag_len);
		for (i = 0; i < m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.%d", filename, i);
			if ((wfd[i] = open(tmpname, O_RDONLY)) < 0) {
				DBG("open('%s') error, skip it....", tmpname);
				ebi->frag_err_list[ebi->nerrs++] = i;
				continue;
			}
			//dbg("ebi->frag_ptrs[%d], len:[%d]", i, ebi->frag_len);
			if(readn(wfd[i], ebi->frag_ptrs[i], ebi->frag_len) != ebi->frag_len) {
				ERR_SYS("readn() error");
			}
			close(wfd[i]);
		}
		if (ebi->nerrs > p) {
			release_ec_buf(ebi);
			err_quit("Too many(%d) fragments lost, must be less(or equal) than [%d], quit", ebi->nerrs, p);
		}

		//printf("####################################\n");
		//dbg("total [%d] fragment lost:", ebi->nerrs);
		//for (i = 0; i < ebi->nerrs; i++) {
		//	printf(" %d ", ebi->frag_err_list[i]);
		//}	
		//printf("\n####################################\n");

		if (ebi->nerrs > 0) {
			int	ret;
			gettimeofday(&start, NULL);
			
			gf_gen_cauchy1_matrix(ebi->encode_matrix, ebi->m, ebi->k);
			ret = gf_gen_decode_matrix_simple(ebi->encode_matrix, ebi->decode_matrix,
					 ebi->invert_matrix, ebi->temp_matrix, ebi->decode_index,
                                          ebi->frag_err_list, ebi->nerrs, k, m);
			if (ret != 0) {
				ERR_QUIT("Fail on generate decode matrix, quit, ret[%d]", ret);
			}
			//for (i = 0; i < k; i++) {
			//	dbg("decode_index[%d]: %d", i, ebi->decode_index[i]);
			//}
			// Pack recovery array pointers as list of valid fragments
			for (i = 0; i < k; i++) {
				ebi->recover_srcs[i] = ebi->frag_ptrs[ebi->decode_index[i]];
			}
			// Recover data
			//ec_init_tables(k, p, &encode_matrix[k * k], g_tbls);
			ec_init_tables(k, ebi->nerrs, ebi->decode_matrix, ebi->g_tbls);
			//ec_encode_data(len, k, nerrs, g_tbls, recover_srcs, recover_outp);
			ec_encode_data(ebi->frag_len, ebi->k, ebi->nerrs, ebi->g_tbls, ebi->recover_srcs, ebi->recover_outp);
		}

		// save frag buffer , don't memcpy()
		for (i = 0; i < ebi->m; i++) {
			ebi->recover_frag_ptrs[i] = ebi->frag_ptrs[i];
		}
		for (i = 0; i < ebi->nerrs; i++) {
			ebi->recover_frag_ptrs[ebi->frag_err_list[i]] = ebi->recover_outp[i];
		}

		msg("####### Recovery time: %ld (us) #########", time_since(&start));
		snprintf(tmpname, sizeof(tmpname), "%s", filename);
		if ((tmpfd = open(tmpname, O_CREAT | O_RDWR, 0644)) < 0) {
			ERR_SYS("open('%s') error", tmpname);
		}
		for (i = 0; i < ebi->k; i++) {
			if(writen(tmpfd, ebi->recover_frag_ptrs[i], ebi->frag_len) != ebi->frag_len) {
				ERR_SYS("writen() error");
			}
		}
		close(tmpfd);
		release_ec_buf(ebi);
		dbg("decoder ok");
	}
	return 0;
}
