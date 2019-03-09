CC= g++
CFLAGS= -std=c++11 -fopenmp -O3 -march=native

SRCS = $(wildcard *.c)
CPP_SRCS = $(wildcard *.cpp)

PROGS = $(patsubst %.c,%,$(SRCS))
CPP_PROGS = $(patsubst %.cpp,%,$(CPP_SRCS))

all: $(PROGS) $(CPP_PROGS)

%: %.c
	$(CC) $(CFLAGS)  -o $@ $<
%: %.cpp
	$(CC) $(CFLAGS)  -o $@ $<
clean:
	rm -f $(PROGS) $(CPP_PROGS)