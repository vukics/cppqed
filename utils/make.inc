# -*- Makefile -*-

optimization = yes

# CXX = g++-4.2

# The flags recommended in GSL 1.8 Reference Manual 

CPPFLAGS = -I. $(foreach dir,$(INCDIR),-I$(dir)) \
	 -g \
	 -fno-common -DLINUX -DGSL_CBLAS #\
	 # -Wall -Wextra -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings \
	 # -Wnon-virtual-dtor -Wshadow -Wsign-conversion -Wno-ignored-qualifiers 

CC = g++
LDFLAGS = -g $(foreach dir,$(LIBDIR),-L$(dir)) -Wl,-rpath,/usr/local/lib

ifeq ($(optimization),yes)
CPPFLAGS += -O2 -DNDEBUG
LDFLAGS  += -O2 -DNDEBUG
else
CPPFLAGS += -DBZ_DEBUG
endif

# In some systems the addition of -latlas may be needed here

LDLIBS = -lm -lgsl -lgslcblas -lblitz -lflens -llapack -lblas
