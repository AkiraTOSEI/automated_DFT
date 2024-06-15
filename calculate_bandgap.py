import pandas as pd
import argparse

def analyze_dos(e_fermi=None, eps=1E-3):
    if e_fermi is None:
        raise ValueError("Fermi energy must be specified.")
    
    # Load data
    dos = pd.read_csv('tmp_DOSCAR', delim_whitespace=True, header=None)
    dos = dos.rename(columns={0: 'Energy', 1: 'DOS(UP)', 2: 'DOS(DOWN)', 3: 'Integ_DOS(UP)', 4: 'Integ_DOS(DOWN)'})
    
    # Calculate if energy values are above the Fermi level
    dos['AboveEf'] = dos['Energy'] - float(e_fermi) > 0
    
    # Get the energy, DOS(UP), DOS(DOWN) of the first point above the Fermi level
    energy1, dup1, ddown1 = dos[dos['AboveEf']].head(4)[['Energy', 'DOS(UP)', 'DOS(DOWN)']].values.tolist()[0]
    energy2, dup2, ddown2 = dos[dos['AboveEf']].head(4)[['Energy', 'DOS(UP)', 'DOS(DOWN)']].values.tolist()[1]
    energy3, dup3, ddown3 = dos[dos['AboveEf']].head(4)[['Energy', 'DOS(UP)', 'DOS(DOWN)']].values.tolist()[2]
    energy4, dup4, ddown4 = dos[dos['AboveEf']].head(4)[['Energy', 'DOS(UP)', 'DOS(DOWN)']].values.tolist()[3]

    if (dup1 > eps and dup2 > eps and dup3>eps and dup4>0) or (ddown1 > eps and ddown2 > eps and ddown3 > eps and ddown4 > 0):
        # どうあがいても金属
        non_metal = False
        band_gap = 0
        print('/home/afujii/coding/automated_DFT/poscar_list/mp-1007907.vasp' ,f', {band_gap:.3f}, {non_metal}')
        return
    elif (0 < dup1 <= eps) and (dup2 ==0) and (dup3 == 0) and (dup4 == 0) and (0 < ddown1 < eps) and (ddown2 ==0) and (ddown3 == 0) and (ddown4 == 0):
        # まぁギリ許容範囲
        dos['AboveEf'] = dos['Energy'] - energy1 > 0
        non_metal = 'True-G'       
    elif (0 < dup1 <= eps) and (0 < dup2 <= eps) and (dup3 == 0) and (dup4 == 0) and (0 < ddown1 < eps) and (0 < ddown2 < eps) and (ddown3 == 0) and (ddown4 == 0):
        # ギリギリ許容範囲
        dos['AboveEf'] = dos['Energy'] - energy2 > 0
        non_metal = 'True-GG'       
    elif (0 < dup1 <= eps) and (0 < dup2 <= eps) and (0 < dup3 <= eps) and (dup4 == 0) and (0 < ddown1 < eps) and (0 < ddown2 < eps) and (0 < ddown3 < eps) and (ddown4 == 0):
        # まじのギリギリのギリ許容範囲
        dos['AboveEf'] = dos['Energy'] - energy4 > 0
        non_metal = 'True-GGG' 
    else:
        dos['AboveEf'] = dos['Energy'] - energy3 > 0
        non_metal = True
    # dup1 < eps and ddown < eps の場合は許容する
    
    # Get Energy of the top of the valence band
    top_VB_idx = dos[dos['AboveEf']].index[0] - 1
    top_VB_energy = dos['Energy'].loc[top_VB_idx]
    
    # Get Energy of the bottom of the conduction band
    dos_above_Ef = dos[dos['AboveEf']][['Energy', 'DOS(UP)', 'DOS(DOWN)']]
    dos_above_Ef['state_exist'] = (dos_above_Ef[['DOS(UP)', 'DOS(DOWN)']] > 0).any(axis=1)
    bottom_CB_energy = dos_above_Ef[dos_above_Ef['state_exist']].head(1)['Energy'].values[0]
    
    band_gap = bottom_CB_energy - top_VB_energy
 
    print('/home/afujii/coding/automated_DFT/poscar_list/mp-1007907.vasp' ,f', {band_gap:.3f}, {non_metal}')

def main():
    parser = argparse.ArgumentParser(description="Analyze DOS data for band gaps.")
    parser.add_argument("-f", "--fermi", type=float, required=True, help="Fermi energy level")
    
    args = parser.parse_args()
    analyze_dos(args.fermi)

if __name__ == "__main__":
    main()
