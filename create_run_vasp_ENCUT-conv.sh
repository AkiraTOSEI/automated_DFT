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





#################################
### ENCUT Convergence ##########
#################################

# 収束判定を初期化
ENCUT_result=FALSE

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
echo "$cutoff \$energy" > encut_test.dat
echo "\$energy" > encut_e.dat
echo "$cutoff" >> encut_cutoff.dat
EOF
    else
        cat >> run_vasp_ENCUT-conv.sh <<EOF
if [ "\$ENCUT_result" == "True" ]; then
    echo "ENCUT Converged."
else
    echo "$cutoff \$energy" >> encut_test.dat
    echo "\$energy" >> encut_e.dat
    echo "$cutoff" >> encut_cutoff.dat
    ENCUT_result=\$(./convergence_checker_EN encut_e.dat)
fi
EOF
    fi
    i=$((i+1))

done

cat >> run_vasp_ENCUT-conv.sh <<EOF
hostname
date
EOF
echo "run_vasp_ENCUT-conv.sh file has been successfully created"

