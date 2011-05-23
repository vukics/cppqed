SHELL = /bin/sh

static = no

SOURCES := $(wildcard *.cc)

EXECS := $(SOURCES:.cc=)

INCDIR := ../include ../include/details
LIBDIR := ../lib /usr/lib/atlas

include ../make.inc

LDLIBS += -lutils

ifeq ($(static),yes)
LDFLAGS += -static
endif

all : $(EXECS)
.PHONY : all clean

# ************
# Main Targets
# ************

# ... are deduced.

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
	@echo "Removing backup files..."
	rm -f *~
	@echo "Removing object files..."
	rm -f *.o
	@echo "Removing dependency files..."
	rm -f *.d*
	@echo "Removing executables..."
	rm -f $(EXECS)
