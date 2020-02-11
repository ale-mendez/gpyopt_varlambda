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
      echo "initer: $ii, maxevals: $jj"
      sed -i "s/^initer: .*$/initer: $ii/g" inpvar_gpyopt.yml
      sed -i "s/^maxevals: .*$/maxevals: $jj/g" inpvar_gpyopt.yml
      python var_gpyopt.py
   done
done

