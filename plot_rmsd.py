import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    holo = mda.Universe('ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.prmtop', 'ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep_prod.nc')
    apo = mda.Universe('../ABL1_cidi_AF_superposed_apo_prep.prmtop', '../ABL1_cidi_AF_superposed_apo_prep_prod.nc')
    ref = mda.Universe('ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.prmtop', 'ABL1_cidi_AF_superposed_2f4j_lig_xtalwat_prep.inpcrd')

    sns.set_theme()
    R_apo = rms.RMSD(apo,
                 ref,
                 select='protein',
                 ref_frame=0)
    R_apo.run()
    df_apo = pd.DataFrame(R_apo.rmsd,
                  columns=['Frame', 'Time (ns)',
                           'RMSD'])

    R_holo = rms.RMSD(holo,
                 ref,
                 select='protein',
                 ref_frame=0)
    R_holo.run()
    df_holo = pd.DataFrame(R_holo.rmsd,
                  columns=['Frame', 'Time (ns)',
                           'RMSD'])

    df_holo['System'] = ['ABL1 AF2 Holo']*len(df_holo)
    df_apo['System'] = ['ABL1 AF2 Apo']*len(df_apo) 
    df = pd.concat((df_holo, df_apo) )
    df['Time (ns)'] = df[ 'Time (ns)']/1000
    sns.lineplot(data=df, x='Time (ns)', y='RMSD', hue='System')
    plt.title('ABL1 AF2 MD simulations RMSD')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD ($\AA$)')

    plt.savefig('holo_apo_rmsd.png')

if __name__ == '__main__':
    main()
