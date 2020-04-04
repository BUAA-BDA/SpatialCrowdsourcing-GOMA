##
#	Makefile for `dna` application
##
CC = gcc
CXX = g++
LIBS = -lrt
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3
TARGET = ./all

all: OPT Ext-GRT Greedy TGOA TGOA-OP TGOA-Greedy 
Ext-GRT: monitor.h input.h Ext-GRT.cpp
	$(CXX) $(CFLAGS) -o Ext-GRT Ext-GRT.cpp $(LIBS) -D WATCH_MEM

OPT: monitor.h input.h OPT.cpp
	$(CXX) $(CFLAGS) -o OPT OPT.cpp $(LIBS) -D WATCH_MEM

Greedy: monitor.h input.h Greedy.cpp
	$(CXX) $(CFLAGS) -o Greedy Greedy.cpp $(LIBS) -D WATCH_MEM

TGOA: monitor.h input.h TGOA.cpp
	$(CXX) $(CFLAGS) -o TGOA TGOA.cpp $(LIBS) -D WATCH_MEM

TGOA-OP: monitor.h input.h TGOA-OP.cpp
	$(CXX) $(CFLAGS) -o TGOA-OP TGOA-OP.cpp $(LIBS) -D WATCH_MEM

TGOA-Greedy: monitor.h input.h TGOA-Greedy.cpp
	$(CXX) $(CFLAGS) -o TGOA-Greedy TGOA-Greedy.cpp $(LIBS) -D WATCH_MEM

.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm OPT Ext-GRT Greedy TGOA TGOA-OP TGOA-Greedy 2> /dev/null || true

