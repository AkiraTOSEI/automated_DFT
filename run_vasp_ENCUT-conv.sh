#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N POSCAR k-convergence
cd $PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=$(cat $PBS_NODEFILE | wc -l)




echo "########################################################"
echo "########################################################"
echo "###  VASP calculation for ENCUT convergence starts! ###"
echo "########################################################"
echo "########################################################"

cp POSCAR_distorted POSCAR

# 収束判定を初期化

cp POSCAR_distorted POSCAR

### create KPOINTS file
echo "k-point initial"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
tail -n 1 kpoints_candidate.csv | tr ',' ' ' >> KPOINTS
echo "0 0 0" >> KPOINTS

### create INCAR file
# INCAR fo ENCUT: 300
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 300" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "$energy" > energy_EN.dat
echo "300" > ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR > force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' > pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 350
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 350" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "350 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "350" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 400
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 400" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "400 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "400" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 450
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 450" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "450 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "450" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 500
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 500" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "500 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "500" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 550
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 550" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "550 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "550" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 600
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 600" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "600 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "600" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 650
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 650" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "650 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "650" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 700
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 700" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "700 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "700" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 750
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 750" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "750 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "750" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 800
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 800" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "800 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "800" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 850
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 850" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "850 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "850" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 900
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 900" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "900 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "900" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 950
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 950" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "950 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "950" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat
### create INCAR file
# INCAR fo ENCUT: 1000
echo "ENCUT POSCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = 1000" >> INCAR
echo "ISMEAR = 0; SIGMA = 0.1" >> INCAR

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free ' OUTCAR | tail -1 | awk '{print $5}')
echo "1000 $energy" >> encut_test.dat
echo "$energy" >> energy_EN.dat
echo "1000" >> ENCUT.dat
awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR >> force_EN.dat 
grep "external pressure" OUTCAR | awk '{print $4}' >> pressure_EN.dat

###### ENCUTの収束判定を行う.
#収束判定ができるcppファイルをコンパイルする 
g++ -std=c++11 -o check_conv convergence_checker.cpp
g++ -std=c++11 -o calculate_l2_norm l2_norm.cpp
g++ -std=c++11 -o convergence_summary convergence_summary.cpp

# 原子数を取得する
NUM_ATOM=$(sed -n '7p' POSCAR | awk '{sum=0; for (i=1; i<=NF; i++) sum+=$i; print sum}')

# Pressureの最終値との差分が5 GPa, 絶対値が5 GPa以下になるようにする. 
# VASPでは圧力単位がKilo Barなので、GPaにへんかんするために0.1をかける
./check_conv pressure_EN.dat 0.1 0.5 last > pressure_diff.dat
./check_conv pressure_EN.dat 0.1 5. zero > pressure_abs.dat
# エネルギー/原子の最終値との差分が1meV以下になるようにする
./check_conv energy_EN.dat "1000/$NUM_ATOM" 1 last > energy_EN_diff.dat
# 力はl2_normをとり、その値が
./calculate_l2_norm force_EN.dat > force_norm.dat
./check_conv force_norm.dat 1. 0.01 last > force_diff.dat

# 結果ファイルを作成する
paste force_diff.dat pressure_abs.dat pressure_diff.dat energy_EN_diff.dat > __tmp_conv.dat
echo "Fx Fy Fz F_norm(eV/Å) F_diff(eV/Å) pressure(GPa) P_diff(GPa) Energy(meV/atom) E_diff(meV/atom)" > __encut_history.dat
awk '{print $1, $2, $7,$8, $10, $11}' __tmp_conv.dat > tmp.dat; paste force_EN.dat tmp.dat| expand | tr -s "[:space:]" >> __encut_history.dat
echo  "ENCUT" > en_tmp.dat; cat ENCUT.dat >> en_tmp.dat 
paste en_tmp.dat __encut_history.dat | expand | tr -s "[:space:]" >  encut_history.dat
# 収束判定のファイルを作成する
awk '{print $3, $6, $9, $12}' __tmp_conv.dat > convergence_EN.dat
BEST_ENCUT=$(./convergence_summary ENCUT.dat convergence_EN.dat )
echo  "ENCUT = $BEST_ENCUT" > BEST_ENCUT.dat
rm __*.dat

# BEST_ENCUTファイルの存在、内容、及び数値のチェック
if [ ! -f BEST_ENCUT.dat ] || [ ! -s BEST_ENCUT.dat ] || ! grep -q '[0-9]' BEST_ENCUT.dat; then
  echo "エラー: 'BEST_ENCUT.dat' ファイルが存在しないか、中身が空、または数値が含まれていません。"
  exit 1
else
  echo "'BEST_ENCUT.dat' ファイルは正常に存在し、中身に数値が含まれています。"
  # ここにBEST_ENCUTファイル用の処理を書く
fi

hostname
date


##### 収束した値に従い、INCARを修正する
# INCARファイルの内容を定義
echo "INCAR with BEST ENCUT" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = $BEST_ENCUT" >> INCAR
echo "ISMEAR = -5; SIGMA = 0.1" >> INCAR
echo "PREC = accurate" >> INCAR
echo "EDIFF = 1.0e-6 # default:10^-4, SCF計算の収束条件,推奨は1E-6らしい。" >> INCAR
echo "ISYM = -1" >> INCAR 
echo "SYMPREC = 0.00000001" >> INCAR

echo "INCAR file has been successfully created."

echo "#######################################################"
echo "###  VASP calculation for ENCUT convergence ends!   ###"
echo "#######################################################"
