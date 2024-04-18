# generate_kpoints.py
import numpy as np
import pandas as pd
import sys

def reciprocal_lattice_vectors(file_path: str) -> np.ndarray:
    with open(file_path, 'r') as file:
        lines = file.readlines()
        a_vector = np.array(list(map(float, lines[2].split())))
        b_vector = np.array(list(map(float, lines[3].split())))
        c_vector = np.array(list(map(float, lines[4].split())))

    a_length = np.linalg.norm(a_vector)
    b_length = np.linalg.norm(b_vector)
    c_length = np.linalg.norm(c_vector)

    reciprocal_lengths = np.array([1/a_length, 1/b_length, 1/c_length])
    return reciprocal_lengths

def calculate_kpoints_adjusted(file_path: str) -> list:
    reciprocal_vectors = reciprocal_lattice_vectors(file_path)
    kpoints_list = []
    for multiplier in range(1, 200):
        multiplied_vectors = reciprocal_vectors * multiplier
        if all(value > 1 for value in multiplied_vectors):
            kpoints_list.append(multiplied_vectors.tolist())
    return kpoints_list


def list_element_check(a:list,b:list):
    assert len(a) == len(b)
    same_elem_exist = False
    for i in range(len(a)):
        if a[i] == b[i]:
            same_elem_exist = True
    return same_elem_exist


def calculate_rounded_kpoints(file_path: str) -> list:
    kpoints_list = calculate_kpoints_adjusted(file_path)
    rounded_kpoints_list = [np.round(kpoint).astype(int).tolist() for kpoint in kpoints_list]
    # 通常の収束判定用のファイル
    new_list = []
    for kpoints in rounded_kpoints_list:
        if not all(value >= 1 for value in kpoints):
            continue
        if len(new_list) == 0:
            new_list.append(kpoints)
        elif not kpoints == new_list[-1]:
            #elif not list_element_check(kpoints, new_list[-1]):
            new_list.append(kpoints)
    
    
    # 収束しなかった時用に、正方晶のKPOINTSを用意する
    new_list2 = []
    for kpoints in rounded_kpoints_list:
        k = kpoints[0]
        if not all(value >= 1 for value in kpoints):
            continue
        if len(new_list2) == 0:
            new_list2.append([k, k, k])
        elif not [k, k, k] == new_list2[-1]:
            #elif not list_element_check(kpoints, new_list[-1]):
            new_list2.append([k, k, k])
    
    
    return new_list[:12]



def main():
    if len(sys.argv) < 1:
        print("Usage: python generate_kpoints.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    rounded_kpoints_list = calculate_rounded_kpoints(file_path)
    
    rounded_kpoints_df = pd.DataFrame(rounded_kpoints_list)
    rounded_kpoints_df.to_csv("kpoints_candidate.csv", index=False,header=False)
    print("K-points candidates have been saved to kpoints_candidate.csv")

if __name__ == "__main__":
    main()

