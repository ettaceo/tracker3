CC = $(CCPREFIX)gcc
AR = $(CCPREFIX)ar
CFLAGS=-Wall -Werror -Wno-unused-function -std=gnu99 -O2

all : algebra3 bisector scenario solution

%.o : %.c
	$(CC) $(CFLAGS) -fpic -o $@ -c $<

algebra3 : algebra3.c algebra3.h
	$(CC) -DLOCAL_BUILD $(CFLAGS) -o $@ $^ -lm

bisector : bisector.o algebra3.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

scenario : scenario.o algebra3.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

solution : solution.o algebra3.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f algebra3 bisector scenario solution *.o
