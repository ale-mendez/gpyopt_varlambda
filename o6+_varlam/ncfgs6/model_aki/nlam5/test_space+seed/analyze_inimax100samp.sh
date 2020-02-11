#!/bin/bash

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
      fmain="i"$ii"_m"$jj
      echo "initer: $ii, maxevals: $jj"
      cd $fmain
      if [ ! -d semillas ]; then	
         mkdir semillas
         mv seed* semillas/.
      fi
      ../analyze_data.py
      cd ../
   done
done

if [ -f var_gpyopt.py ]; then
   rm autovarlambda* NIST_* var_gpyopt.py 
fi
