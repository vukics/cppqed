SHELL = /bin/sh

SRCDIRS := quantumdata quantumoperator structure quantumtrajectory elements/frees elements/interactions elements/composites utils/src

EXECS := $(patsubst scripts/%.cc,%,$(wildcard scripts/*.cc))

# Scripts to be put into directory "scripts" and are then found automatically

vpath %.cc scripts:$(SRCDIRS)
vpath %.h $(SRCDIRS):utils/include

include make.macosx.inc

optimization = yes
profiling = no

CXX = g++

CPPFLAGS = -ftemplate-depth-128 -w -fPIC $(foreach dir,$(INCDIR),-I$(dir)) -DGSL_CBLAS

LDFLAGSLIBS := $(foreach dir,$(LIBDIR),-L$(dir))

ifeq ($(optimization),yes)
CPPFLAGS += -finline-functions    -O3 -DNDEBUG
else
CPPFLAGS += -fno-inline        -g -O0 -DBZ_DEBUG
endif

LDLIBSTHIRDPARTY := -lgsl -lgslcblas -lblitz
LDLIBS = $(LDLIBSTHIRDPARTY) -lC++QED

ifeq ($(with-flens),no)
CPPFLAGS += -DDO_NOT_USE_FLENS
EXCLUDED := utils/src/DrivenDampedHarmonicOscillator.o # If this is included in libC++QED.so, then linkage without FLENS will not work
else
LDLIBS += -lflens 
endif

ifeq ($(profiling),yes)
CPPFLAGS += -pg
LDFLAGS += -pg
endif


all : $(EXECS)
.PHONY : all clean library


# ************
# Main Targets
# ************

STDOBJ := $(filter-out $(EXCLUDED),$(foreach dir,$(SRCDIRS),$(patsubst $(dir)/%.cc,$(dir)/%.o,$(wildcard $(dir)/*.cc))))

$(EXECS) : $(CPPQEDLIB)

$(CPPQEDLIB) : $(STDOBJ)
	g++ $(LDFLAGSLIBS) $(LDLIBSTHIRDPARTY) $(CPPQEDLIBLINKFLAGS) -o $(CPPQEDLIB) $(STDOBJ)

clean :
	@echo "Removing library file..."
	rm -f $(CPPQEDLIB)
	@echo "Removing object files..."
	rm -f $(STDOBJ)
	@echo "Removing executables..."
	rm -f $(EXECS)
