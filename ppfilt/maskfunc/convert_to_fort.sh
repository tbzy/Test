
sed 's/E/d/g' maskr.15 | awk '{ print "    mfunc(", NR-1, ") = ", $2 }'





