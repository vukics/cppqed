rm $1.0*.jpg
latex $1; dvips -E -i -S 1 $1
# for the pdf version:
# for file in $1.0*; do convert -density 1200 $file $file.jpg; rm $file; done

for file in $1.0*; do convert -density 180 $file $file.jpg; rm $file; done