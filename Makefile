SHELL = /bin/sh

static = no

SRCDIRS := quantumdata quantumoperator quantumtrajectory structure elements/frees elements/interactions elements/composites

EXECS := $(patsubst scripts/%.cc,%,$(wildcard scripts/*.cc))

# Driver programs directory "driver" and are then found automatically

EXCLUDED := MovingParticleCavity.o ParticleAlongCavity.o MovingParticleBase.o
# These are excluded from compiling & linking of an EXEC
# Nevertheless, they can still be compiled by invoking explicitely make foo.o

vpath %.cc scripts:$(SRCDIRS)
vpath %.h $(SRCDIRS):C++Utils/include C++Utils/include/details

INCDIR := $(SRCDIRS) C++Utils/include C++Utils/include/details
LIBDIR := C++Utils/lib /usr/lib/atlas

include C++Utils/make.inc

LDLIBS += -lC++Utils

ifeq ($(static),yes)
LDFLAGS += -static
endif

LIBUTILS := C++Utils/lib/libC++Utils.a

all : $(EXECS)
.PHONY : all utils install clean

install : all ; mv $(EXECS) ~/bin/

$(LIBUTILS) :
	$(MAKE) -CC++Utils s

utils : $(LIBUTILS)

# ************
# Main Targets
# ************
STDOBJ := $(filter-out $(EXCLUDED),$(foreach dir,$(SRCDIRS),$(patsubst $(dir)/%.cc,%.o,$(wildcard $(dir)/*.cc))))

$(EXECS) : $(STDOBJ) $(LIBUTILS)

# ************************
# Determining dependencies 
# ************************
%.d: %.cc
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

SOURCES := $(foreach dir,$(SRCDIRS) scripts,$(wildcard $(dir)/*.cc))
-include $(SOURCES:.cc=.d)

clean :
	@echo "Removing backup files..."
	rm -f *~ scripts/*~ $(foreach dir,$(SRCDIRS),$(dir)/*~)
	@echo "Removing object files..."
	rm -f *.o
	@echo "Removing dependency files..."
	rm -f scripts/*.d* test/*.d* $(foreach dir,$(SRCDIRS),$(dir)/*.d*)
	@echo "Removing precompiled headers..."
	rm -f $(foreach dir,$(SRCDIRS),$(dir)/*.h.gch)
	@echo "Removing executables..."
	rm -f $(EXECS)
	make -C C++Utils clean
