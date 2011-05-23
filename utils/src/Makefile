SHELL = /bin/sh
INCDIR := ../include .
LIBDIR := ../lib

EXCLUDED := Combinatorics.o

ARCH := $(LIBDIR)/libC++Utils.a

include ../make.inc

SOURCES := $(wildcard *.cc)
STDOBJ := $(filter-out $(EXCLUDED),$(patsubst %.cc,%.o,$(SOURCES)))

all : $(ARCH)($(STDOBJ))
.PHONY : all clean

# ************************
# Determining dependencies 
# ************************

%.d: %.cc
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(SOURCES:.cc=.d)

clean :
	rm -f *~
	rm -f *.d*
	rm -f *.o
	rm -f $(ARCH)
