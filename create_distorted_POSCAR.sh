#!/bin/bash

cat > distort_poscar.py <<EOF
import numpy as np

# POSCARファイルのパスを指定
file_path = 'POSCAR_gt'  # このパスを実際のファイルパスに置き換えてください

# ファイルを読み込む
with open(file_path, 'r') as file:
    lines = file.readlines()

# 2行目から4行目の結晶ベクトルを読み込む
crystal_vectors = np.array([line.split() for line in lines[2:5]], dtype=float)
num_total_atoms = np.array(lines[6].strip().split(" ")).astype(int).sum()
atom_vectors = np.array([line.split()[:3] for line in lines[8:] if line.strip()], dtype=float)
assert len(atom_vectors) == num_total_atoms

total_distortion = 0.1 # Angstrom

random_distortion_size_for_vec0 = (1.-np.random.uniform(size=num_total_atoms))*total_distortion
norm_vec0 = np.sqrt(np.sum(np.square(crystal_vectors[0])))
random_distortion_size_for_vec1 = (total_distortion-random_distortion_size_for_vec0)*np.random.uniform(size=num_total_atoms)
norm_vec1 = np.sqrt(np.sum(np.square(crystal_vectors[1])))
random_distortion_size_for_vec2 = total_distortion-random_distortion_size_for_vec1-random_distortion_size_for_vec0
norm_vec2 = np.sqrt(np.sum(np.square(crystal_vectors[2])))
direct_distortion_vectors = np.stack([
    random_distortion_size_for_vec0 / norm_vec0, 
    random_distortion_size_for_vec1 / norm_vec1, 
    random_distortion_size_for_vec2 / norm_vec2, 
]).T

(np.sqrt(np.sum(np.square(np.matmul(direct_distortion_vectors, crystal_vectors)),axis=1)) < total_distortion).all()


# POSCARファイルの1行目から7行目をそのまま保持し、8行目以降をatom_vectorsを使用して再構成するためのコードを作成します。

# 既に読み込んでいるlines変数から1行目から7行目を抽出
header_lines = lines[:8]

# atom_vectorsを文字列に変換してファイルに書き込むための形式にします
atom_vectors_str = '\n'.join([' '.join(map(str, vector)) for vector in (atom_vectors+direct_distortion_vectors)])

# 新しいPOSCAR内容を組み立てます
new_poscar_content = ''.join(header_lines)  + atom_vectors_str

# 新しいPOSCARファイルを書き出します
new_poscar_path = 'POSCAR_distorted'
with open(new_poscar_path, 'w') as new_file:
    new_file.write(new_poscar_content)

EOF

cp POSCAR POSCAR_gt

python distort_poscar.py

cp POSCAR_distorted POSCAR
