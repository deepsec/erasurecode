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

#include "common.h"
#include "error.h"

#define PAGE_SIZE		4096
#define DEFAULT_IO_SIZE		(PAGE_SIZE * 1024)		// 4M
#define MAX_IO_SIZE		(PAGE_SIZE * 1024 * 128)	// 256M
#define MIN_IO_SIZE		(PAGE_SIZE * 128)		// 512K
#define MAX_IOCB_COUNT		8192
#define SHARDING_SIZE		(1024 * 1024 * 1024 * 4L);	// 4G

typedef struct thread_sharding_info
{
	int	infd;
	int	outfd;
	loff_t	offset;
	size_t	size;
	int64_t	io_size;
} SHARDING_INFO;


void libaio_read_prepare(int fd, struct iocb **iocb_list, int iocb_list_len, loff_t io_size, loff_t offset, io_callback_t cb)
{
	struct iocb *iocb_ptr;
	void *buf;
	loff_t i;

	for (i = 0; i < iocb_list_len; ++i) {
		if ((iocb_ptr = malloc(sizeof(struct iocb))) == NULL) {
			ERR_SYS("malloc(%d) error", sizeof(struct iocb));
		}
		if (posix_memalign(&buf, PAGE_SIZE, io_size) != 0) {
			ERR_QUIT("posix_memalign(align_size='%d', io_size='%d') error", PAGE_SIZE, io_size);
		}
		io_prep_pread(iocb_ptr, fd, buf, io_size, offset+i*io_size);
		io_set_callback(iocb_ptr, cb);
		iocb_list[i] = iocb_ptr;
	}
}

void file_sharding_aio_copy(void *arg)
{
	SHARDING_INFO	*si = (SHARDING_INFO *)arg;
	char		*buf;
	ssize_t		n;
	int		rfd, wfd;
	uint64_t	io_blocks, last;
	uint64_t	curpos, io_size;

	struct io_event	*r_events, *r_event_ptr, *w_events, *w_event_ptr;
	struct iocb	**r_iocb_list, **w_iocb_list;
	struct iocb 	*rp, *wp;
	uint64_t	iocb_list_len;
	io_context_t	r_ctx, w_ctx;

	rfd = si->infd;
	wfd = si->outfd;
	io_size = si->io_size;

	io_blocks = si->size / io_size;
	last = si->size % io_size;	
	iocb_list_len = io_blocks;	/* prepare libaio iocbs */

	if (iocb_list_len > MAX_IOCB_COUNT) {
		ERR_QUIT("iocb_list_len[%d] is too big, quit", iocb_list_len);
	}
	//dbg("start copy: offset[%ld], size[%ld], io_size[%ld] ......", si->offset, si->size, io_size);
	
	if (iocb_list_len > 0) {	// libaio engine for PAGE_SIZE aligned data
		int	alldone, done, j;

		set_fl(rfd, O_DIRECT);
		set_fl(wfd, O_DIRECT);
		curpos = si->offset;

		memset(&r_ctx, 0, sizeof(r_ctx));
		if (io_setup(iocb_list_len, &r_ctx) != 0) {
			ERR_SYS("io_setup('%d') error", iocb_list_len);
		}
		if ((r_iocb_list = malloc(sizeof(struct iocb *) * iocb_list_len)) == NULL) {
			ERR_SYS("malloc(%d) error", sizeof(struct iocb *));
		}
		libaio_read_prepare(rfd, r_iocb_list, iocb_list_len, io_size, curpos, NULL);
		if (io_submit(r_ctx, iocb_list_len, r_iocb_list) != iocb_list_len) {
			ERR_SYS("io_submit() error");
		}
		r_events = malloc(sizeof(struct io_event) * iocb_list_len);
		if (r_events == NULL) {
			ERR_SYS("malloc() error");
		}

		memset(&w_ctx, 0, sizeof(w_ctx));
		if (io_setup(iocb_list_len, &w_ctx) != 0) {
			ERR_SYS("io_setup('%d') error", iocb_list_len);
		}
		if ((w_iocb_list = malloc(sizeof(struct iocb *) * iocb_list_len)) == NULL) {
			ERR_SYS("malloc(%d) error", sizeof(struct iocb *));
		}
		w_events = malloc(sizeof(struct io_event) * iocb_list_len);
		if (w_events == NULL) {
			ERR_SYS("malloc() error");
		}

		alldone = 0;
		while (1) {
			if ((done = io_getevents(r_ctx, 1, iocb_list_len - alldone, r_events + alldone, NULL)) < 1) {
				ERR_SYS("io_getevents() error");
			}
			for (j=0, r_event_ptr=r_events+alldone; j < done; ++j) {
				if (r_event_ptr[j].res2 != 0) {
					ERR_SYS("aio pread error()");
				}
				rp = r_event_ptr[j].obj;
				wp = malloc(sizeof(struct iocb));
				if (wp == NULL) {
					ERR_SYS("malloc() error");
				}
				io_prep_pwrite(wp, wfd, rp->u.c.buf, rp->u.c.nbytes, rp->u.c.offset);
				io_set_callback(wp, NULL);
				w_iocb_list[j+alldone] = wp;
			}
			if (io_submit(w_ctx, done, w_iocb_list + alldone) != done) {
				ERR_SYS("io_submit() error");
			}
			//dbg("********************************* alldone:[%d]  done: [%d] ******************************", alldone, done);
			alldone += done;
			if (alldone == iocb_list_len) {
				break;
			}
		}

		alldone = 0;
		while (1) {
			if ((done = io_getevents(w_ctx, 1, iocb_list_len - alldone, w_events + alldone, NULL)) < 1) {
				ERR_SYS("io_getevents() error");
			}
			for (j=0, w_event_ptr=w_events+alldone; j < done; ++j) {
				if (w_event_ptr[j].res2 != 0) {
					ERR_SYS("aio pwrite error()");
				}
			}
			alldone += done;
			if (alldone == iocb_list_len) {
				break;
			}
		}	
		for (j = 0; j < iocb_list_len; j++) {
			io_callback_t cb  = (io_callback_t)w_events[j].data;
			wp = w_events[j].obj;
			if (cb) {
				cb(w_ctx, wp, w_events[j].res, w_events[j].res2);
			}
			free(wp->u.c.buf);
		}
		for (j = 0; j < iocb_list_len; j++) {
			free(r_iocb_list[j]);
			free(w_iocb_list[j]);
		}
		free(r_iocb_list);
		free(w_iocb_list);
		free(r_events);
		free(w_events);
		io_destroy(r_ctx);
		io_destroy(w_ctx);
	}

	if (last > 0) { /* the last no-memaligned must be use Buffer IO */
		dbg("enter the last no-memaligned copy");
		buf = malloc(io_size);
		if (NULL == buf) {
			ERR_SYS("malloc(%d) error", io_size);
		}
		curpos = si->offset + io_size * iocb_list_len;
		clear_fl(rfd, O_DIRECT);
		clear_fl(wfd, O_DIRECT);
		n = preadn(rfd, buf, io_size, curpos);
		if (n < 0) {
			ERR_SYS("preadn() error");
		}
		if (pwriten(wfd, buf, n, curpos) < 0) {
			ERR_SYS("pwriten() error");
		}
		free(buf);
	}
	//dbg("finish copy: offset[%ld], size[%ld], io_size[%ld] ......", si->offset, si->size, io_size);
}

void copy_file_attribute(char *src, char *dst)
{
	struct stat	s_st, d_st;
	uid_t		uid;
	gid_t		gid;
	mode_t		mode;
	struct utimbuf	t;
	
	if (stat(src, &s_st) < 0 || lstat(dst, &d_st)) {
		ERR_SYS("stat('%s') or lstat('%s') error", src, dst);
	}
	t.actime = s_st.st_atime;
	t.modtime = s_st.st_mtime;
	mode = s_st.st_mode;
	uid = s_st.st_uid;
	gid = s_st.st_gid;
	
	/* ignore error */
	utime(dst, &t);
	chmod(dst, mode);
	chown(dst, uid, gid);

	return;
}

int set_symlink_timestamp(const char *pathname, struct timespec atime, struct timespec mtime)
{
	int	ret;
	struct timespec	ts[2];
	ts[0] = atime;
	ts[1] = mtime;
	ret = utimensat(AT_FDCWD, pathname, ts, AT_SYMLINK_NOFOLLOW);
	return ret;
}
	
int main(int argc, char **argv)
{
	int infd, outfd;
	size_t file_size, file_shardings, file_lastsharding;
	struct stat	src_st, dst_st;
	SHARDING_INFO	*sinfo;
	int		sinfo_counts;
	int64_t		i, total_copy;
	int64_t		cmd_io_size, io_size, sharding_size;
	int		opt;
	int		keep_attr;
	char		src[NAME_MAX], dst[NAME_MAX];
	struct timeval	tv_begin, tv_end;
	time_t		time_elapsed;

	cmd_io_size = DEFAULT_IO_SIZE / 1024;
	keep_attr = 0;
	while ((opt = getopt(argc, argv, "i:k")) != -1)
        {
                switch (opt)
                {
                        case 'i':
				cmd_io_size = strtoul(optarg, NULL, 10);
				break;
                        case 'k':
                               	keep_attr = 1;	
                                break;
                        default:
				err_quit("USAGE: %s [-i io_size(KB)] [-k] <src_file> <dst_file>", argv[0]);
                }
        }

	/* transfer into KB(for io_size) and  MB(for sharding sharding_size) */
	io_size = cmd_io_size * 1024;
	sharding_size = SHARDING_SIZE;
	if (io_size > MAX_IO_SIZE) {
		dbg("cmdline io_size[%d KB] is great than [%d KB], set it to [%d KB]", cmd_io_size, MAX_IO_SIZE/1024, MAX_IO_SIZE/1024);
		io_size = MAX_IO_SIZE;
	} else if (io_size < MIN_IO_SIZE) {
		dbg("cmdline io_size[%d KB] is less than [%d KB], set it to [%d KB]", cmd_io_size, MIN_IO_SIZE/1024, MIN_IO_SIZE/1024);
		io_size = MIN_IO_SIZE;
	}

	if (argc - optind != 2) {
		err_quit("USAGE: %s [-i io_size(KB)] [-k] <src_file> <dst_file>", argv[0]);
	}
	strncpy(src, argv[optind], sizeof(src));
	strncpy(dst, argv[optind+1], sizeof(dst));

	if ((infd = open(src, O_RDONLY)) < 0) {
		ERR_SYS("open() %s error", argv[1]);
	}
	if (lstat(src, &src_st) < 0) {
		ERR_SYS("lstat('%s') error", src);
	}
	/* only deal with regular file and symlink */
	if (!(S_ISREG(src_st.st_mode) || S_ISLNK(src_st.st_mode))) {
		ERR_QUIT("src file[%s] isn't regular file, also not symlink, quit", src);
	}

	if ((outfd = open(dst, O_CREAT | O_WRONLY, 0644)) < 0) {
		if (EISDIR == errno) {
                	strncat(dst, "/", sizeof(dst));
                	strncat(dst, src, sizeof(dst));
			if ((outfd = open(dst, O_CREAT | O_WRONLY, 0644)) < 0) {
				ERR_SYS("open() %s error", dst);
			}
        	}
		else {
			ERR_SYS("open() %s error", dst);
		}
	}
	if (lstat(dst, &dst_st) < 0) {
		ERR_SYS("lstat('%s') error", dst);
	}
	/* if same file, quit */
	if (dst_st.st_ino == src_st.st_ino) {
		ERR_QUIT("src['%s'], dst['%s'] is same file", src, dst);
	}
	/* if symlink, copy symlink, and return */
	if (S_ISLNK(src_st.st_mode)) {
		char	linkname[NAME_MAX];

		unlink(dst);
		if (readlink(src, linkname, sizeof(linkname)) < 0) {
			ERR_SYS("readlink() error");
		}
		if (symlink(linkname, dst) < 0) {
			ERR_SYS("symlink() error");
		}
		if (1 == keep_attr) {
			set_symlink_timestamp(dst, src_st.st_atim, src_st.st_mtim);
		}
		return 0;
	}

	if (ftruncate(outfd, 0) < 0) {
		ERR_SYS("ftruncate() error");
	}

	file_size = src_st.st_size;
	file_shardings = file_size / sharding_size;
	file_lastsharding = file_size % sharding_size;

	dbg("src_file[%s], file_size:[%lld] io_size[%lld], sharding_size[%lld], dst_file[%s]", src, file_size, io_size, sharding_size, dst);

	if (file_lastsharding > 0) { //the last is not zero
		sinfo_counts = file_shardings + 1;
	} else {
		sinfo_counts = file_shardings;
	}

	sinfo = malloc(sinfo_counts * sizeof(SHARDING_INFO));
	if (sinfo == NULL) {
		ERR_SYS("sinfo malloc() error");
	}
	memset(sinfo, 0, sinfo_counts * sizeof(SHARDING_INFO));
	dbg("file_shardings:[%d], sinfo_counts[%d]", file_shardings, sinfo_counts);
	for (i = 0; i < sinfo_counts; i++) {
		sinfo[i].infd = infd;
		sinfo[i].outfd = outfd;
		sinfo[i].offset = sharding_size * i;
		sinfo[i].size = sharding_size;
		sinfo[i].io_size = io_size;
		if ((i == (sinfo_counts - 1)) && (file_lastsharding > 0)) {	// the last sharding less than the sharding_size 
				sinfo[i].size = file_lastsharding;
		}
	}

	if (gettimeofday(&tv_begin, NULL) < 0) {
		ERR_SYS("gettimeofday() error");
	}
	total_copy = 0;
	for (i = 0; i < sinfo_counts; i++) {
		file_sharding_aio_copy(&sinfo[i]);
		total_copy += sinfo[i].size;
		dbg("sharding: [%d], copied: [%.2f MB] ---  [%2.2f%%]", i, (total_copy*1.0)/(1024*1024), total_copy*100.0/file_size);
	}	
	if (gettimeofday(&tv_end, NULL) < 0) {
		ERR_SYS("gettimeofday() error");
	}
	time_elapsed = tv_end.tv_sec - tv_begin.tv_sec;
	if (time_elapsed == 0) {
		time_elapsed = 1;
	}
	dbg("Total copied: [%lld bytes], elapsed: [%lld secs],  Speed: [%.2f MB/s]", total_copy, time_elapsed, (total_copy*1.0)/(1024.0*1024.0)/(time_elapsed*1.0));
	free(sinfo);
	close(infd);
	close(outfd);

	if (1 == keep_attr) {
		copy_file_attribute(src, dst);
	}
	dbg("*****  COPY DONE *****");

	return 0;
}
