#!/bin/bash

~/autovarlambda/gpyopt_varlambda/o6+_varlam/link_all.sh

init=$(date)
echo $init >> info_run.out
for i in {1..1}
do
   loop_init=$(date)
   echo $i $loop_init >> info_run.out
   fseed="seed_"$i 
   iseed=$(( $i*1234+($i-1)*56789 ))
   echo "seed_"$i 
   mkdir $fseed
   sed -i "s/^iseed: .*$/iseed: $iseed/g" inpvar_gpyopt.yml
   python ./var_gpyopt.py 
   mv O* $fseed/.
done
fin=$(date)
echo $fin >> info_run.out

rm autovarlambda* NIST_* var_gpyopt.py 
