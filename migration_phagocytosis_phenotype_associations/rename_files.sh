#!/bin/bash
for NAME in $(cut -f1 ../../data/migration_assay_11_n_12_transwells.txt);do

echo $NAME
FIXED_NAME=$(echo $NAME | grep -o '[^/]*$')
echo $FIXED_NAME
NEW_NAME=$(cat ../../data/migration_assay_11_n_12_transwells.txt | grep $FIXED_NAME | awk '{print $2}')
echo $NEW_NAME
mv ../../data/crams/$FIXED_NAME.cram ../../data/crams/$NEW_NAME.cram
done
