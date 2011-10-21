CPPQEDLIB := libC++QED.dylib
CPPQEDLIBLINKFLAGS = -dynamiclib -install_name `pwd`/libC++QED.dylib
INCDIR := $(SRCDIRS) utils/include /opt/local/include
LIBDIR := /opt/local/lib
with-flens = no
LDFLAGS = -L. $(LDFLAGSLIBS)
