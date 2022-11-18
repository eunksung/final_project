#!/bin/bash
# read any tsv files in the directory you are working on
for i in *.tsv; do
# print the file without first 7 lines and take 4th column of the file
tail -n +7 "$i" | awk '{print $4}' > "$i.tmp" && mv -f "$i.tmp" "$i"
# give a name to the column that we have taken out with its file name 
( echo -e "${i%.*}"; cat $i; ) > ${i%.*}.txt
done
