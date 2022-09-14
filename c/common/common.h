#ifndef __COMMON_H__
#define __COMMON_H__

ssize_t readn(int fd, void *ptr, size_t n);
ssize_t writen(int fd, const void *ptr, size_t n);
ssize_t preadn(int fd, void *ptr, size_t n, loff_t offset);
ssize_t pwriten(int fd, const void *ptr, size_t n, loff_t offset);
int	get_fl(int fd);
void	set_fl(int fd, int flags);
void	clear_fl(int fd, int flags);

#endif
