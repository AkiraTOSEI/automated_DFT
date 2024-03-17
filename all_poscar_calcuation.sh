#!/bin/bash

# poscar_dataディレクトリ内の各ファイルに対してループ
for file in poscar_list/*; do
    # all_convergence.shスクリプトを実行し、現在のファイル名を引数として渡す
    sh all_convergence.sh "$file"
done

