#!/bin/bash

init=$(date)
echo $init >> info_run.out
for i in {1..100}
do
   loop_init=$(date)
   echo $i $loop_init >> info_run.out
   fseed="seed_"$i 
   iseed=$(( $i*1234+($i-1)*56789 ))
   echo "seed_"$i 
   mkdir $fseed
   sed -i "s/^iseed: .*$/iseed: $iseed/g" inpvar_gpyopt.yml
#   tail -n1 inpvar_gpyopt.yml > $fseed/seed.inp
   ./varpol_gpyopt.py 
   mv Be* $fseed/.
done
rm adas* adf* oic ols olg OMG* TERMS LEVELS CONFIG*
fin=$(date)
echo $fin >> info_run.out
