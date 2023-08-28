#!/bin/bash

for i in */; do
    if [ -f $i/descriptors.json ]; then
        if [ -d $i/modeller_veryfast_0 ]; then
            cp $i/modeller_veryfast_0/sequence.seq $i               
            cp $i/modeller_veryfast_0/sequence.seq.ali $i
            cp $i/modeller_veryfast_0/complex.pdb $i/complex0.pdb
            rm -fr $i/modeller_veryfast_*
        fi
    fi
done
