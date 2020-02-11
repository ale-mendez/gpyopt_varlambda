#!/bin/bash

~/autovarlambda/gpyopt_varlambda/varpol/link_all.sh
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
   python ./varpol_gpyopt.py 
   mv Be* $fseed/.
done
fin=$(date)
echo $fin >> info_run.out
./xtractdata.sh
rm autovarlambda* das* NIST_* varpol_gpyopt.py xtractdata.sh
