import numpy as np

# POSCARファイルのパスを指定
file_path = 'POSCAR_gt'  # このパスを実際のファイルパスに置き換えてください

# ファイルを読み込む
with open(file_path, 'r') as file:
    lines = file.readlines()

kpoint_density = 4000

# 2行目から4行目の結晶ベクトルを読み込む
crystal_vectors = np.array([line.split() for line in lines[2:5]], dtype=float)
num_total_atoms = np.array(lines[6].strip().split(" ")).astype(int).sum()
atom_vectors = np.array([line.split()[:3] for line in lines[8:] if line.strip()], dtype=float)
assert len(atom_vectors) == num_total_atoms     

# K点の密度に沿って、KPOINTSファイルに書き込むファイルを作る。
crystal_vectors_length = np.sqrt(np.sum(np.square(crystal_vectors),axis=1))
norm_ratio = (1./crystal_vectors_length) /np.max(crystal_vectors_length )
base_val = np.power(kpoint_density/num_total_atoms/np.prod(norm_ratio), 1./3.)
best_kpoints = np.round(base_val*norm_ratio).astype(np.int64)

# ISMEAR=-5 を用いるため、KPOINTSの値は4より大きくなくてはならない。もしくは、再計算されたK点密度が小さすぎてもいけない
recalculated_k_density = np.prod(best_kpoints) * num_total_atoms
if recalculated_k_density < 4000*0.75 or np.min(best_kpoints)<=3:
    k_strings = ''
else:
    k_strings = ''
    for k in best_kpoints:
        k_strings += f'{k} '
    k_strings = k_strings.strip()
print(k_strings)
