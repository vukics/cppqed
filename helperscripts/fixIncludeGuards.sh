# Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#!/bin/bash

for f in $(find * -name '*.h' -o -name '*.tcc'); do perl -pi -e "s[(\S*_INCLUDED)]{my \$newguard=uc(qq(${f}_INCLUDED)); \$newguard=~ s|[./]|_|g; qq(\$newguard)}e" $f; done
