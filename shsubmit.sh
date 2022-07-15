#!/bin/bash
for x in $(seq 0 20)
do
y=$((x*28+1))

var=`bash submit.sh $y`

echo $var

done