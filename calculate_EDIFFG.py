import numpy as np

# POSCARファイルのパスを指定
file_path = 'POSCAR_gt'  # このパスを実際のファイルパスに置き換えてください

# ファイルを読み込む
with open(file_path, 'r') as file:
    lines = file.readlines()

    base_EDIFF = -1E-3

# 2行目から4行目の結晶ベクトルを読み込む
crystal_vectors = np.array([line.split() for line in lines[2:5]], dtype=float)
num_total_atoms = np.array(lines[6].strip().split(" ")).astype(int).sum()
atom_vectors = np.array([line.split() for line in lines[8:] if line.strip()], dtype=float)
assert len(atom_vectors) == num_total_atoms    

EDIFFG = base_EDIFFG * num_total_atoms

print(f'EDIFFG = {EDIFFG:.2e}')
