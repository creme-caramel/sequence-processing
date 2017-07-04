FILEPATH := $(realpath $(lastword $(MAKEFILE_LIST)))
CURDIR := $(shell cd $(dir $(FILEPATH));pwd)
SRCDIR := $(CURDIR)/src/

CXX = g++
CFLAGS = -Wall -g -std=c++11
BIN = main

all: $(BIN)

$(BIN):
	$(CXX) $(CFLAGS) -o $(BIN) $(SRCDIR)main.cpp $(SRCDIR)groupedseqs_a.cpp

clean:
	rm -f $(BIN)

.PHONY:
	clean all
