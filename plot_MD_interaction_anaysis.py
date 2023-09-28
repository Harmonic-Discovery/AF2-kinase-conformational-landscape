import MDAnalysis as mda
import prolif as plf
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# This code was created using examples from prolif

def main():
    # Load trajectory
    system = sys.argv[1]
    u = mda.Universe(system+'.prmtop', system+'_prod.nc')
    lig = u.select_atoms("resname VX6")
    prot = u.select_atoms("protein")

    #Analyze interactions 
    fp = plf.Fingerprint()
    fp.run(u.trajectory, lig, prot)
    df = fp.to_dataframe()
    df = df.droplevel("ligand", axis=1)
    occ = df.mean()
    occ_df = occ.reset_index()
    occ_df = occ_df.rename(columns={0: 'frequency'})

    # Create the bar plot using Seaborn
    plt.figure(figsize=(10, 6))
    occ_df['residue number'] = [x[:3]+str(int(x[3:]) + 240) for x in occ_df['protein'].values]
    sns.barplot(x='residue number', y='frequency', hue='interaction', data=occ_df)
    plt.xticks(rotation=45)
    plt.title('Interaction analysis of ABL1 AF2 structure with 2F4J ligand during MD simulation')
    plt.legend(bbox_to_anchor=(1.2, 1), borderaxespad=0)
    plt.tight_layout()
    plt.savefig(system+'_interactions.png')


if __name__ == '__main__':
    main()
