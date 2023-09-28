import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def prepare_df(system, num_actives, where, to_ignore=None, best_poses=None):
    """
    Picks only the best pose for each ligand and label the ligands.
    The actives are required to come first in the docked file, followed by inactives/decoys
    """
    mols = open(where+'/'+system+'_smina-docked.sdf', "r") 
    print(' Loaded docking output file')
    lines = mols.readlines()
    str_to_check = 'minimizedAffinity'

    check = 0
    affinities = []
    for line in lines:
        ls = line.split()
        if check == 1:
            affinities.append(ls[0])

        if 'min' in line:
            check = 1
        else:
            check = 0

    if to_ignore is None:
        actives = pd.DataFrame(affinities[:num_actives], columns=['minimizedAffinity'])

    else:
        aff_act_np  = np.array(affinities[:num_actives*9]).reshape(num_actives, 9)# smina generated 9 poses
        aff_act_np[to_ignore] = 0
        aff_act_best = [aff_act_np[x,y] for x,y in zip(range(len(aff_act_np[:,0])), best_poses-1)]
        actives = pd.DataFrame(aff_act_best, columns=['minimizedAffinity'])

    inactives_decoys = pd.DataFrame(affinities[num_actives::], columns=['minimizedAffinity'])
    actives['type'] = 'active'
    inactives_decoys['type'] = 'inactives+decoys'

    df = pd.concat((actives, inactives_decoys))
    return df


def calculate_enrichments(df):
    """
    Sorts the dataframe using the docking binding affinity and calculate the enrichment
    """

    df['minimizedAffinity'] = df['minimizedAffinity'].astype(float)
    df_sorted = df.sort_values(by=['minimizedAffinity'])
    
    efs_act = [] 
    perc_tot = []
    tot_act = len(df_sorted[df_sorted['type'] == 'active'])
    tot_comp = len(df_sorted)

    for i in range(1, len(df_sorted)+1):
        df_sub = df_sorted[:i]
        num_act = len(df_sub[df_sub['type'] == 'active'])

        ef_act = num_act / tot_act * 100
        perc = i/tot_comp * 100
       
        efs_act.append(ef_act)
        perc_tot.append(perc)
        
    return efs_act, perc_tot


def plot_efs_var(target, systems, top, apo, af_sys, states, ef_act_tot, ef_dec_tot, log, df_auc):
    """
    Plot enrichment curves in linear or semilogarithmic scale
    """
    for i in range(len(states)):
        ef_act = ef_act_tot[i]
        ef_dec = ef_dec_tot[i]
        y_min = [1000]*len(ef_act[0])
        y_max = [-10]*len(ef_act[0])
        state_systems = systems
        plt.figure()
        plt.tight_layout()
        sns.set_context('notebook', font_scale=1)
        print('Number of systems in this state: '+str(len(state_systems)))
        for m in range(len(state_systems)):
            y_sys = ef_act[m]
            auc = df_auc[df_auc['System'] == state_systems[m]]['AUC']
            print('number of data point: '+str(len(y_sys)))
            
            if state_systems[m] != af_sys:
                for x in range(len(y_sys)):
                    if y_sys[x] <= y_min[x]:
                        y_min[x] = y_sys[x]
                    if y_sys[x] >= y_max[x]:
                        y_max[x] = y_sys[x]

            if state_systems[m] == top:
                sns.lineplot(x=ef_dec[m], y=ef_act[m], label=state_systems[m]+'_holo, LogAUC = '+str(round(auc.values[0],2)))
            elif state_systems[m] == apo:
                sns.lineplot(x=ef_dec[m], y=ef_act[m], label=state_systems[m]+'_apo, LogAUC = '+str(round(auc.values[0],2)))
            elif state_systems[m] == af_sys:
                sns.lineplot(x=ef_dec[m], y=ef_act[m], label=state_systems[m]+', LogAUC = '+str(round(auc.values[0], 2)))
        
        plt.fill_between(ef_dec[0], y_min, y_max, alpha=0.2)

        liny = np.linspace(0,100, len(ef_dec[m]))
        linx = np.linspace(0,100, len(ef_act[m]))
        sns.lineplot(x=linx, y=liny, linestyle='dashed', color='gray')

        plt.title('Docking Enrichment for '+target+' '+states[i].upper()+' states')
        plt.ylabel('% of active compounds')
        plt.axvspan(0, 10, alpha=.3, color='gray')

        if log == True:
            plt.xscale('log')
            plt.xlabel('%  of total compounds (log)')
            plt.tight_layout()
            plt.savefig(target+'_'+states[i]+'_log-enrichment.jpg', dpi=300)

        else:
            plt.xlabel('%  of total compounds')
            plt.tight_layout()
            plt.savefig(target+'_'+states[i]+'-enrichment.jpg', dpi=300)

        plt.close()

