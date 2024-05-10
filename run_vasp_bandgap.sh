#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd $PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=$(cat $PBS_NODEFILE | wc -l)



#################################
######     bandgap       ##########
#################################


# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p bandgap_cal
cp POTCAR KPOINTS calculate_bandgap.py ./bandgap_cal
cd bandgap_cal

echo "################################################"
echo "################################################"
echo "###  VASP calculation for  Bandgap starts!  ###"
echo "################################################"
echo "################################################"
echo ""

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_KPOINTS=$(cat ../BEST_KPOINTS.dat)
echo "BEST KPOINTS"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "$BEST_KPOINTS" >> KPOINTS
echo "0 0 0" >> KPOINTS


# ENCUTの収束後値をもとに、bandgap計算用のINCARを作成する.
BEST_ENCUT=$(cat ../BEST_ENCUT.dat)
echo "bandap calculation INCAR" > INCAR
echo "ISTART = 1 #default:1, 0:scratch, 1:WAVECARを読む、なければ0になる " >> INCAR
echo "ICHARG = 11 #default:2(ISTART=0)0,else:0, 0:初期波動関数から電荷を計算,2:原子電荷の重ね合わせを使う,11:CHGCAR利用" >> INCAR
echo "ISPIN = 2 # default:1, 1:スピン分極なし、2:スピン分極あり" >> INCAR
echo "$BEST_ENCUT" >> INCAR
echo "EDIFF = 1.0e-6 # default:10^-4, SCF計算の収束条件,推奨は1E-6らしい。" >> INCAR
echo "LWAVE = .F. # WAVECARを出力するか" >> INCAR
echo "ISMEAR = -5 # default:1,部分占有軌道の設定、絶縁体・半導体では-5を推奨 " >> INCAR
echo "#SIGMA = 0.01 # default:0.2, smearingの幅、系の事前知識がないときはISMEAR=0でSIGMA=0.03~0.05 " >> INCAR
echo "PREC = accurate" >> INCAR
echo "LREAL = .F. # 精度の指定。default NORMAL" >> INCAR
echo "NELM = 200 # default:60, SC計算の最大回数を指定" >> INCAR
#echo "NPAR = 4 # defalut:コア数、並列計算に使うコア数" >> INCAR
echo "LORBIT = 11 # default: None, 11:DOSCAR, PROCARを出力" >> INCAR
echo "NEDOS = 2000 #default:301, DOS計算のグリッド数" >> INCAR
cat ../INCAR_tail >> INCAR

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4

echo "#######################################"
echo "###            VASP ends!           ###"
echo "#######################################"

cp OUTCAR ../final_OUTCAR
cp DOSCAR ../final_DOSCAR
cp POSCAR ../final_POSCAR
cp KPOINTS ../final_KPOINTS
cp INCAR ../final_INCAR

### バンドギャップ計算用のDOSCARを作成する
# 最初の8行をスキップし、4列以上のデータを含む行を無視
tail -n +9 DOSCAR | awk 'NF <= 5 { print }' > tmp_DOSCAR

E_fermi=$(grep "E-fermi" ../final_OUTCAR | awk '{print $3}')
OUTPUT_FILE="../../bandgap_result.csv"

echo "E_fermi: $E_fermi"

# outputが存在するかどうかをチェック
if [ -e $OUTPUT_FILE ]; then
    # ファイルが存在する場合、通常のlsの結果をファイルに出力
    python calculate_bandgap.py -f $E_fermi >> $OUTPUT_FILE
else
    # ファイルが存在しない場合、隠しファイルを含むls -aの結果をファイルに出力
    echo 'file_name, bandgap, nonmetal' > $OUTPUT_FILE
    python calculate_bandgap.py -f $E_fermi >> $OUTPUT_FILE
fi

echo "#######################################"
echo "###     bandgap cpp ends!           ###"
echo "#######################################"

cd ../

# 特定のファイル以外を削除
#rm -rf bandgap_cal final_sc relax_cal
#shopt -s extglob
#rm -f !(POSCAR_gt|POSCAR_distorted|final_INCAR|final_KPOINTS|final_OUTCAR|final_DOSCAR|final_POSCAR|used_POTCAR.txt|kp_history.dat|encut_history.dat)

date



