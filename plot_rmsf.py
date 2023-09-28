import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-sys', '--system', type=str,  required=True, help='')

    args = parser.parse_args()
    system = args.system
    sns.set_theme()

    holo = mda.Universe('ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.prmtop', 'ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep_prod.nc')
    apo = mda.Universe('../ABL1_cidi_AF_superposed_apo_prep.prmtop', '../ABL1_cidi_AF_superposed_apo_prep_prod.nc')
    ref = mda.Universe('ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.prmtop', 'ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.inpcrd')

    # Compute rmsf
    protein = 'protein'
    aligner = align.AlignTraj(holo, ref,
                          select='not resname ACE and not resname NME and protein and name CA',
                          in_memory=True).run()
    c_alphas = holo.select_atoms('not resname ACE and not resname NME and protein and name CA')
    R_holo = rms.RMSF(c_alphas).run()
    
    aligner = align.AlignTraj(apo, ref,
                          select='not resname ACE and not resname NME and protein and name CA',
                          in_memory=True).run()
    c_alphas = apo.select_atoms('not resname ACE and not resname NME and protein and name CA')
    R_apo = rms.RMSF(c_alphas).run()

    # Adjust the residue numbering
    resids =  c_alphas.resids+240
    df = pd.DataFrame()
    df['RMSF']  = np.concatenate((R_holo.results.rmsf,R_apo.results.rmsf))
    apo_label = ['ABL1 AF2 Apo']*len(R_apo.results.rmsf)
    holo_label = ['ABL1 AF2 Holo']*len(R_holo.results.rmsf)
    df['System'] = holo_label + apo_label
    df['Resnum'] = np.concatenate((resids, resids))
    df.to_pickle('ABL1_cidi_AF_RMSF.pkl')

    # plot
    plt.title('ABL1 AF2 MD simulations backbone RMSF')
    sns.lineplot(data=df, x='Resnum', y='RMSF', hue='System')
    plt.xlabel('Residue number')
    plt.ylabel('RMSF ($\AA$)')
    plt.savefig(system+'_rmsf.png')

    # make pdb with rmsf at b-factor for rendering
    holo.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
    protein = holo.select_atoms('protein') # select protein atoms
    for residue, r_value in zip(protein.residues, R_holo.rmsf):
        residue.atoms.tempfactors = r_value
    holo.atoms.write('rmsf_tempfactors.pdb')

if __name__ == '__main__':
    main()
