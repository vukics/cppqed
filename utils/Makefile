SHELL = /bin/sh

SUBS = s t ts

all : $(SUBS)
.PHONY : all s t ts clean pack

s :
	$(MAKE) -C src

t :
	$(MAKE) -C test

ts :
	$(MAKE) -C testsuite

clean :
	$(MAKE) clean -C src
	$(MAKE) clean -C test
	$(MAKE) clean -C testsuite
	rm -f *~ include/*~ include/*.h.gch

pack : clean
	tar -C.. -czf ../C++Utils.tgz C++Utils

