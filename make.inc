CPPQEDLIB := libC++QED.so
CPPQEDLIBLINKFLAGS := -shared
INCDIR := $(SRCDIRS) utils/include
LIBDIR :=
with-flens = yes
LDFLAGS = -L. $(LDFLAGSLIBS) -Wl,-R -Wl,"`pwd`" -Wl,-rpath-link -Wl,"`pwd`"
