#!/bin/bash

if [ -f total_loss.dat ]; then
    rm total_loss.dat
fi
if [ -f params.dat ]; then
    rm params.dat
fi

for i in {1..100}
do
    grep -i "Total loss" seed_$i/Be*.out >> total_loss.dat
    grep -i "alpha" seed_$i/Be*.out >> params.dat
    grep -i "rcut" seed_$i/Be*.out >> params.dat
done

