#!/bin/bash

# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <VASP_FILE_PATH>"
    exit 1
fi
FILE_REL_PATH="$1"
PRECISION="${2:-high}" # precision引数がなければ、デフォルトで'high'を使用
KP_CONV_THRES="${3:-0.03}" # precision引数がなければ、デフォルトで'0.05'を使用
EN_CONV_THRES="${3:-0.0001}" # precision引数がなければ、デフォルトで'0.05'を使用

# ファイルの存在を確認
if [ ! -f "$FILE_REL_PATH" ]; then
    echo "Error: File does not exist - $FILE_REL_PATH"
    exit 2
fi

# パスを絶対パスに変換する
FILE_PATH=$(realpath "$FILE_REL_PATH")

# working dirctoryを作ってそこで作業を行う
WORK_DIR="workspace__"
WORK_DIR="$WORK_DIR""$(basename "$FILE_PATH")"
mkdir -p $WORK_DIR
cp create_conv_cpp.sh create_INCAR.sh create_POSCAR.sh create_POTCAR.sh generate_kpoints.py ./$WORK_DIR
cp create_run_vasp_KP-conv.sh create_run_vasp_ENCUT-conv.sh create_run_vasp_relax.sh sh create_run_vasp_aft_relax_sc.sh create_run_vasp_bandgap_cal.sh calculate_bandgap.sh ./$WORK_DIR
cd $WORK_DIR

# 収束判定のファイルを作る
sh create_conv_cpp.sh $KP_CONV_THRES
g++ -o convergence_checker_KP convergence_checker.cpp
sh create_conv_cpp.sh $EN_CONV_THRES
g++ -o convergence_checker_EN convergence_checker.cpp


sh create_INCAR.sh
sh create_POSCAR.sh $FILE_PATH
sh create_POTCAR.sh $FILE_PATH $PRECISION
# K点の候補ファイルを作る
python generate_kpoints.py  $FILE_PATH 

### scriptを書く
# K点の候補ファイルをもとに、全部を計算するスクリプトを書く
sh create_run_vasp_KP-conv.sh  $FILE_PATH
# ENCUTの収束のスクリプトを書く
sh create_run_vasp_ENCUT-conv.sh $FILE_PATH
# 緩和計算のスクリプトを書く
sh create_run_vasp_relax.sh
# 緩和計算後のsc計算のスクリプトを書く
sh create_run_vasp_aft_relax_sc.sh
# bandgap計算のスクリプトを書く(OUTCARからの抽出も含む)
sh create_run_vasp_bandgap_cal.sh $FILE_PATH


# vaspの計算scriptをつなげる
cat run_vasp_KP-conv.sh > run_vasp_all.sh
tail -n +8 run_vasp_ENCUT-conv.sh >> run_vasp_all.sh
tail -n +8 run_vasp_relax.sh >> run_vasp_all.sh
tail -n +8 run_vasp_aft_relax_sc.sh >> run_vasp_all.sh
tail -n +8 run_vasp_bandgap.sh >> run_vasp_all.sh

# vaspの計算を行う
qsub run_vasp_all.sh
