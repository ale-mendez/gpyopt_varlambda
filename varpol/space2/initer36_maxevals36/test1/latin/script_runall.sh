#!/bin/bash

for i in {1..10}
do
   frun="tmp_run"$i
   mkdir $frun
   cd $frun
   ln -s ~/autovarlambda/gpyopt_varpol/autovarlambda.cpython-37m-x86_64-linux-gnu.so
   ln -s ~/autovarlambda/gpyopt_varpol/das_6CFG
   ln -s ~/autovarlambda/gpyopt_varpol/NIST_cfgs.dat
   ln -s ~/autovarlambda/gpyopt_varpol/NIST_energies.dat
   ln -s ~/autovarlambda/gpyopt_varpol/NIST_lines.dat
   ln -s ~/autovarlambda/gpyopt_varpol/NIST_terms.dat
   ln -s ~/autovarlambda/gpyopt_varpol/varpol_gpyopt.py
   cp ../inpvar_gpyopt.yml .
   cp ~/autovarlambda/gpyopt_varpol/run_gpyopt.batch .
   ini=$(( ($i-1)*10+1 ))
   fin=$(( $i*10 ))
   echo $ini $fin
   sed -i "s/^for i in .*$/for i in {$ini..$fin}/g" run_gpyopt.batch
   qsub -N "i24e24" run_gpyopt.batch
   cd ..
done


