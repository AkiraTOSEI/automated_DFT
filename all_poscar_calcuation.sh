#!/bin/bash

# bandgapの計算ファイルを初期化
echo "file_name, bandgap, nonmetal" > bandgap_result.csv 

# poscar_dataディレクトリ内の各ファイルに対してループ
for file in poscar_list/*; do
    # all_convergence.shスクリプトを実行し、現在のファイル名を引数として渡す
    sh all_convergence.sh "$file"
done

