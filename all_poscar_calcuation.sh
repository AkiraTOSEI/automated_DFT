#!/bin/bash


# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <POSCAR_DIRECTORY_PATH>"
    exit 1
fi
FILE_REL_PATH="$1"

# bandgapの計算ファイルを初期化
echo "file_name, bandgap, nonmetal" > bandgap_result.csv 

# poscar_dataディレクトリ内の各ファイルに対してループ
echo "num poscar files: $(ls $FILE_REL_PATH | wc -l)"
for file in $FILE_REL_PATH/*; do
    # all_convergence.shスクリプトを実行し、現在のファイル名を引数として渡す
    sh all_convergence.sh "$file"
done

