#!/bin/bash

# Check if a file name is provided as the first argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <FILE_NAME>"
  exit 1
fi

FILE_NAME="$1"

cat > run_vasp_ENCUT-conv.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N $FILE_NAME k-convergence
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)


echo "#######################################################"
echo "#######################################################"
echo "###  VASP calculation for ENCUT convergence starts! ###"
echo "#######################################################"
echo "#######################################################"
echo ""

EOF


i=1
for cutoff in 200 250 300 350 400 450 500 550 600 650 700 750 800; do
    cat >> run_vasp_ENCUT-conv.sh <<EOF
if [ "\$ENCUT_result" == "True" ]; then
    echo "ENCUT Converged."
else
    ### create KPOINTS file
    # INCAR fo ENCUT: $cutoff
    echo "ENCUT $FILE_NAME" > INCAR
    echo "ISTART = 0 ; ICHARG = 2" >> INCAR
    echo "ENCUT = $cutoff" >> INCAR
    echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR
   
    ### calculcation
    mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
    energy=\$(grep 'free ' OUTCAR | tail -1 | awk '{print \$5}')
fi
EOF

    if [ $i -eq 1 ]; then
        cat >> run_vasp_ENCUT-conv.sh <<EOF
echo "\$energy" > energy_EN.dat
echo "$cutoff" > ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print \$4, \$5, \$6}' OUTCAR > force.dat 
grep "external pressure" OUTCAR | awk '{print \$4}' > pressure.dat
EOF
    else
        cat >> run_vasp_ENCUT-conv.sh <<EOF
echo "$cutoff \$energy" >> encut_test.dat
echo "\$energy" >> energy_EN.dat
echo "$cutoff" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print \$4, \$5, \$6}' OUTCAR >> force.dat 
grep "external pressure" OUTCAR | awk '{print \$4}' >> pressure.dat
EOF
    fi
    i=$((i+1))
done



cat >> run_vasp_ENCUT-conv.sh <<EOF

###### ENCUTの収束判定を行う.
#収束判定ができるcppファイルをコンパイルする 
g++ -std=c++11 -o check_conv convergence_checker.cpp
g++ -std=c++11 -o calculate_l2_norm l2_norm.cpp
g++ -std=c++11 -o convergence_summary convergence_summary.cpp

# 原子数を取得する
total_lines=\$(wc -l < POSCAR)
NUM_ATOM=\$((total_lines - 8))

# Pressureの最終値との差分が5 GPa, 絶対値が5 GPa以下になるようにする. 
# VASPでは圧力単位がKilo Barなので、GPaにへんかんするために0.1をかける
./check_conv pressure.dat 0.1 0.5 last > pressure_diff.dat
./check_conv pressure.dat 0.1 5. zero > pressure_abs.dat
# エネルギー/原子の最終値との差分が1meV以下になるようにする
./check_conv energy_EN.dat "1000/\$NUM_ATOM" 1 last > energy_EN_diff.dat
# 力はl2_normをとり、その値が
./calculate_l2_norm force.dat > force_norm.dat
./check_conv force_norm.dat 1. 0.01 last > force_diff.dat

# 結果ファイルを作成する
paste force_diff.dat pressure_abs.dat pressure_diff.dat energy_EN_diff.dat > __tmp_conv.dat
echo "Fx Fy Fz F_norm(eV/Å) F_diff(eV/Å) pressure(GPa) P_diff(GPa) Energy(meV/atom) E_diff(meV/atom)" > __encut_history.dat
awk '{print \$1, \$2, \$7,\$8, \$10, \$11}' __tmp_conv.dat > tmp.dat; paste force.dat tmp.dat| expand | tr -s "[:space:]" >> __encut_history.dat
echo  "ENCUT" > en_tmp.dat; cat ENCUT.dat >> en_tmp.dat 
paste en_tmp.dat __encut_history.dat | expand | tr -s "[:space:]" >  encut_history.dat
# 収束判定のファイルを作成する
awk '{print \$3, \$6, \$9, \$12}' __tmp_conv.dat > convergence.dat
BEST_ENCUT=\$(./convergence_summary ENCUT.dat convergence.dat )
echo  \$BEST_ENCUT
rm __*.dat
hostname
date


echo "#######################################################"
echo "###  VASP calculation for ENCUT convergence ends!   ###"
echo "#######################################################"
EOF
echo "run_vasp_ENCUT-conv.sh file has been successfully created"

