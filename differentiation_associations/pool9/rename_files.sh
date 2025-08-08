#!/bin/bash
for NAME in $(cut -f1 ../../data/run45137_sampleInfo.txt);do  

echo $NAME
FIXED_NAME=$(echo $NAME | grep -o '[^/]*$')
echo $FIXED_NAME
NEW_NAME=$(cat ../../data/run45137_sampleInfo.txt | grep -w $FIXED_NAME | awk '{print $2}')
echo $NEW_NAME
mv ../../data/crams/$FIXED_NAME.cram ../../data/crams/$NEW_NAME.cram
done
