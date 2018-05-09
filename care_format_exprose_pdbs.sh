#!/bin/bash

for i in `seq 1 $1`; do
    mv out.pdb.${i} out_${i}.pdb
done
