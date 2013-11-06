#!/bin/bash

for f in $(find * -name '*.h' -o -name '*.tcc'); do perl -pi -e "s[(\S*_INCLUDED)]{my \$newguard=uc(qq(${f}_INCLUDED)); \$newguard=~ s|[./]|_|g; qq(\$newguard)}e" $f; done
