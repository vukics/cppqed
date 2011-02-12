cd _build/html
for i in `find . -name "*.html"`; do 
  echo $i
  sed -i "
	{s/__PL__/+/g}
	{s/__MI__/-/g}
	{s/__CO__/,/g}
	{s/__RP__/)/g}
	{s/__LP__/(/g}
	{s/C++QED vMilestone 8 documentation/C++QEDv2 Milestone8 documentation/g}
	{s/Python Module Index/Index of Header Files/g}
        {
N
s/Note<\/p>\n<p>Template argument definitions:/Template argument definitions/
}
        "\
	$i; 
done
cd ../..