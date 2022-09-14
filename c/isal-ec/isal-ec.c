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

#include "common.h"
#include "error.h"

#define M_K_P_MAX	255
#define K_DEFAULT	6
#define P_DEFAULT	3

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

	dbg("release all buffer");
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

int main(int argc, char **argv)
{
	int             opt;
	int		fd, tmpfd;
	struct stat	st;
	int64_t		file_size, frag_len;
	int		m, k, p;
	int		i, is_decode;
	char		filename[NAME_MAX], tmpname[NAME_MAX];
	EC_BUF_INFO	*ebi;
	struct timeval	start;

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

	strncpy(filename, argv[optind], sizeof(filename));
	if (is_decode == 0) {
		if (lstat(filename, &st) < 0) {
			ERR_SYS("lstat('%s') error", filename);
		}
		if ((fd = open(filename, O_RDONLY)) < 0) {
			ERR_SYS("open('%s') error", filename);
		}
		file_size = st.st_size;
		frag_len = file_size / k;
		if (file_size % k) {
			frag_len += 1;
			dbg("**** padding frag for the last fragment");
		}
		dbg("file['%s'], file_size[%ld], m[%d], k[%d], p[%d], frag_len[%ld]", filename, file_size, m, k, p, frag_len);

		ebi = alloc_ec_buf(m, k, p, frag_len);
		for (i = 0; i < k; i++) {
			if (readn(fd, ebi->frag_ptrs[i], frag_len) != frag_len) {
				ERR_SYS("readn() error)");
			}
		}

		gettimeofday(&start, NULL);
		gf_gen_cauchy1_matrix(ebi->encode_matrix, m, k);
		// Initialize g_tbls from encode matrix
		ec_init_tables(k, p, &(ebi->encode_matrix)[k * k], ebi->g_tbls);
		// Generate EC parity blocks from sources
		ec_encode_data(frag_len, k, p, ebi->g_tbls, ebi->frag_ptrs, &(ebi->frag_ptrs)[k]);

		dbg("cost time: %lld (us)", time_since(&start));
		for (i = 0; i < ebi->m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.isal.%d", filename, i);
			if ((tmpfd = open(tmpname, O_CREAT | O_RDWR, 0644)) < 0) {
				ERR_SYS("open('%s') error", tmpname);
			}
			dbg("write to file: '%s'", tmpname);
			if(writen(tmpfd, ebi->frag_ptrs[i], ebi->frag_len) != ebi->frag_len) {
				ERR_SYS("writen() error");
			}
			close(tmpfd);
		}
		release_ec_buf(ebi);
		close(fd);
	}
	else {
		// get frag_len
		for (i = 0; i < m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.isal.%d", filename, i);
			if (lstat(tmpname, &st) < 0) {
				err_msg("lstat('%s') error, skip it!", tmpname);
				continue;
			}
			frag_len = st.st_size;
			break;
		}
		ebi = alloc_ec_buf(m, k, p, frag_len);
		dbg("m[%d], k[%d], p[%d], frag_len[%d]", m, k, p, frag_len);
		dbg("m[%d], k[%d], p[%d], frag_len[%d]", ebi->m, ebi->k, ebi->p, ebi->frag_len);
		for (i = 0; i < m; i++) {
			snprintf(tmpname, sizeof(tmpname), "%s.isal.%d", filename, i);
			if ((tmpfd = open(tmpname, O_RDONLY)) < 0) {
				err_msg("open('%s') error", tmpname);
				dbg("skip it");
				ebi->frag_err_list[ebi->nerrs++] = i;
				continue;
			}
			dbg("ebi->frag_ptrs[%d], len:[%d]", i, ebi->frag_len);
			if(readn(tmpfd, ebi->frag_ptrs[i], ebi->frag_len) != ebi->frag_len) {
				ERR_SYS("readn() error");
			}
			close(tmpfd);
		}
		if (ebi->nerrs > p) {
			release_ec_buf(ebi);
			err_quit("Too many(%d) fragments lost, must be less(or equal) than [%d], quit", ebi->nerrs, p);
		}
		printf("####################################\n");
		dbg("total [%d] fragment lost:", ebi->nerrs);
		for (i = 0; i < ebi->nerrs; i++) {
			printf(" %d ", ebi->frag_err_list[i]);
		}	
		printf("\n####################################\n");
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
			for (i = 0; i < k; i++) {
				dbg("decode_index[%d]: %d", i, ebi->decode_index[i]);
			}
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

		dbg("####### Recovery time: %ld #########", time_since(&start));
		snprintf(tmpname, sizeof(tmpname), "%s.isal", filename);
		if ((tmpfd = open(tmpname, O_CREAT | O_RDWR, 0644)) < 0) {
			ERR_SYS("open('%s') error", tmpname);
		}
		for (i = 0; i < ebi->k; i++) {
			if(writen(tmpfd, ebi->frag_ptrs[i], ebi->frag_len) != ebi->frag_len) {
				ERR_SYS("writen() error");
			}
		}
		close(tmpfd);
		release_ec_buf(ebi);
		dbg("decoder ok");
	}
	return 0;
}
