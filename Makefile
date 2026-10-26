#  Makefile for interpBilinear.
#
#  The library is header-only; there is nothing to build except the test.

CC     ?= gcc
CFLAGS ?= -O2 -std=c99 -pedantic -Wall -Wextra
LIBS    = -lm

all: test

testInterpBilinear: testInterpBilinear.c interpBilinear.h
	$(CC) $(CFLAGS) -o $@ testInterpBilinear.c $(LIBS)

test: testInterpBilinear
	./testInterpBilinear

clean:
	rm -f testInterpBilinear

.PHONY: all test clean
