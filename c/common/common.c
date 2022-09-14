#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include "error.h"
#include "common.h"

/* Read "n" bytes from a descriptor  */
ssize_t readn(int fd, void *ptr, size_t n)
{
	size_t nleft;
	ssize_t nread;

	nleft = n;
	while (nleft > 0) {
		if ((nread = read(fd, ptr, nleft)) < 0) {
			if (nleft == n)
				return (-1);	/* error, return -1 */
			else
				break;	/* error, return amount read so far */
		}
		else if (nread == 0) {
			break;	/* EOF */
		}
		nleft -= nread;
		ptr += nread;
	}
	return (n - nleft);	/* return >= 0 */
}

/* Write "n" bytes to a descriptor  */
ssize_t writen(int fd, const void *ptr, size_t n)
{
	size_t nleft;
	ssize_t nwritten;

	nleft = n;
	while (nleft > 0) {
		if ((nwritten = write(fd, ptr, nleft)) < 0) {
			if (nleft == n)
				return (-1);	/* error, return -1 */
			else
				break;	/* error, return amount written so far */
		}
		else if (nwritten == 0) {
			break;
		}
		nleft -= nwritten;
		ptr += nwritten;
	}
	return (n - nleft);	/* return >= 0 */
}


/* Read "n" bytes from a descriptor  */
ssize_t preadn(int fd, void *ptr, size_t n, loff_t offset)
{
	size_t nleft;
	ssize_t nread;
	loff_t	curpos;

	nleft = n;
	curpos = offset;
	while (nleft > 0) {
		if ((nread = pread(fd, ptr, nleft, curpos)) < 0) {
			if (nleft == n)
				return (-1);	/* error, return -1 */
			else
				break;	/* error, return amount read so far */
		}
		else if (nread == 0) {
			break;	/* EOF */
		}
		nleft -= nread;
		ptr += nread;
		curpos += nread;
	}
	return (n - nleft);	/* return >= 0 */
}

/* Write "n" bytes to a descriptor  */
ssize_t pwriten(int fd, const void *ptr, size_t n, loff_t offset)
{
	size_t nleft;
	ssize_t nwritten;
	loff_t	curpos;

	nleft = n;
	curpos = offset;
	while (nleft > 0) {
		if ((nwritten = pwrite(fd, ptr, nleft, curpos)) < 0) {
			if (nleft == n)
				return (-1);	/* error, return -1 */
			else
				break;	/* error, return amount written so far */
		}
		else if (nwritten == 0) {
			break;
		}
		nleft -= nwritten;
		ptr += nwritten;
		curpos += nwritten;
	}
	return (n - nleft);	/* return >= 0 */
}

int get_fl(int fd)
{
	int val;

	if ((val = fcntl(fd, F_GETFL, 0)) < 0) {
		ERR_SYS("fcntl F_GETFL error");
	}
	return val;
}

void set_fl(int fd, int flags)
{
	int val;

	if ((val = fcntl(fd, F_GETFL, 0)) < 0) {
		ERR_SYS("fcntl F_GETFL error");
	}
	val |= flags;
	if (fcntl(fd, F_SETFL, val) < 0) {
		ERR_SYS("fcntl F_SETFL error");
	}
}

void clear_fl(int fd, int flags)
{
	int val;

	if ((val = fcntl(fd, F_GETFL, 0)) < 0) {
		ERR_SYS("fcntl F_GETFL error");
	}
	val &= ~flags;
	if (fcntl(fd, F_SETFL, val) < 0) {
		ERR_SYS("fcntl F_SETFL error");
	}
}
