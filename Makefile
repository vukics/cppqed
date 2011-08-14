SHELL = /bin/sh

SRCDIRS := quantumdata quantumoperator structure quantumtrajectory elements/frees elements/interactions elements/composites utils/src

EXECS := $(patsubst scripts/%.cc,%,$(wildcard scripts/*.cc))

# Scripts to be put into directory "scripts" and are then found automatically

vpath %.cc scripts:$(SRCDIRS)
vpath %.h $(SRCDIRS):utils/include

INCDIR := $(SRCDIRS) utils/include

optimization = yes
with-flens = yes

CXX = g++

CPPFLAGS = -ftemplate-depth-128 -w -fPIC $(foreach dir,$(INCDIR),-I$(dir)) -DGSL_CBLAS

LDFLAGS = -L. -Wl,-R -Wl,"`pwd`" -Wl,-rpath-link -Wl,"`pwd`"

ifeq ($(optimization),yes)
CPPFLAGS += -finline-functions    -O3 -DNDEBUG
else
CPPFLAGS += -fno-inline        -g -O0 -DBZ_DEBUG
endif

LDLIBS = -lgsl -lgslcblas -lblitz -lC++QED

ifeq ($(with-flens),no)
CPPFLAGS += -DDO_NOT_USE_FLENS
EXCLUDED := utils/src/DrivenDampedHarmonicOscillator.o # If this is included in libC++QED.so, then linkage without FLENS will not work
else
LDLIBS += -lflens 
endif


all : $(EXECS)
.PHONY : all clean library


# ************
# Main Targets
# ************

STDOBJ := $(filter-out $(EXCLUDED),$(foreach dir,$(SRCDIRS),$(patsubst $(dir)/%.cc,$(dir)/%.o,$(wildcard $(dir)/*.cc))))

$(EXECS) : libC++QED.so

libC++QED.so : $(STDOBJ)
	g++ -shared -o libC++QED.so $(STDOBJ)

clean :
	@echo "Removing library file..."
	rm -f libC++QED.so
	@echo "Removing object files..."
	rm -rf $(STDOBJ)
	@echo "Removing executables..."
	rm -f $(EXECS)
