#!/bin/bash

~/autovarlambda/gpyopt_varlambda/o6+_varlam/link_all.sh

iini=1
ifin=5
 
for (( i=$iini; i<=$ifin; i++ ))
do
   jini=$(($i-1))
   jfin=$(($i+1))
   ii=$(($i*10))
   for (( j=$jini; j<=$jfin; j++ ))
   do
      jj=$(($j*10))
      sed -i "s/^initer: .*$/initer: $ii/g" inpvar_gpyopt.yml
      sed -i "s/^maxevals: .*$/maxevals: $jj/g" inpvar_gpyopt.yml
      fmain="i"$ii"_m"$jj
      echo "initer: $ii, maxevals: $jj"
      mkdir $fmain
      cp inpvar_gpyopt.yml das* $fmain
      cd $fmain
      ~/autovarlambda/gpyopt_varlambda/o6+_varlam/link_all.sh
      for k in {1..100}
      do
         fseed="seed_"$k
         iseed=$(( $k*1234+($k-1)*56789 ))
         echo "seed "$k 
         mkdir $fseed
         sed -i "s/^iseed: .*$/iseed: $iseed/g" inpvar_gpyopt.yml
         python ./var_gpyopt.py 
         mv O* $fseed/.
      done
      cd ../
   done
done

rm autovarlambda* NIST_* var_gpyopt.py 
