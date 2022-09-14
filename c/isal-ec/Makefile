CC := gcc
CFLAGS := -Wall -g -D__DEEPDBG__ -D_FILE_OFFSET_BITS=64
INCLUDE := -I. -I../common -I/usr/local/include 
LIBS := -L/usr/local/lib -L/usr/lib -L. -L../common -lpthread -lisal

COMMON_SRCS := ../common/error.c ../common/common.c
COMMON_OBJS := $(subst .c,.o, $(COMMON_SRCS))

EXEC := isal-ec

all: $(EXEC)
	
isal-ec: isal-ec.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio
AIOCopy: AIOCopy.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio

ABIOCopy: ABIOCopy.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio

aioEpoll: aioEpoll.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio

aiocp: aiocp.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio

test: test.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -laio

testso: test.o
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) -ldeepcommon

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ 





.PHONY: all clean
clean:
	@rm -f *.o lib*.so lib*.a $(EXEC) $(COMMON_OBJS)


## file dependency
../common/error.o: ../common/error.c ../common/error.h
../common/common.o: ../common/common.c ../common/error.h
