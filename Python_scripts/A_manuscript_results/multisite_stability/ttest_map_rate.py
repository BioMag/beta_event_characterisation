"""
Calculate and plot multi site stability of rate
"""
import pandas as pd
import mne
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from config import (fname, MRI_dir, spacing)
from settings_hmm_beta import (task, group_id)
from scipy.stats import ttest_ind

# Read subjects and the characteristics csv
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])
rate_feature_df = pd.read_csv(fname.feature_csv(feature = 'rate',job_id=group_id, task=task))

# Read subjects and the characteristics csv from site 2
df_subjects_site2 = pd.read_csv("subject_text_files/list_of_subjects_site2.txt", names=["subject"])
rate_feature_df_site2 = pd.read_csv(fname.feature_csv(feature = 'rate',job_id=group_id+'_site2', task=task))

# Read the parcels
labels_parc = mne.read_labels_from_annot(subject = 'fsaverage', parc = 'aparc_sub', subjects_dir=MRI_dir)
lh_parc= [ lab for lab in labels_parc if (lab.hemi == 'lh' or lab.hemi == 'rh') ]

# Read the source space
fname_src = fname.src(hmm_bids_dir=fname.hmm_bids_dir(subject='fsaverage', ses='01'), subject='fsaverage', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Make empty dictionaries for
p_values_s1_s2 = {'low':[],'high':[]}

rate_values_all = {'low':{'site1':np.zeros((len(df_subjects), len(lh_parc))),
                        'site2':np.zeros((len(df_subjects_site2), len(lh_parc)))},
                'high':{'site1':np.zeros((len(df_subjects), len(lh_parc))),
                        'site2':np.zeros((len(df_subjects_site2), len(lh_parc)))}}

# Make columns to take from the df
columns = []
for sens in range(0,len(lh_parc)):
    columns.append('sens' + str(sens+1))

for low_high in ["low", "high"]:
    
    for i,row in df_subjects.iterrows():
        subject = row["subject"]
        subject_df = rate_feature_df.loc[
            rate_feature_df['subject']==subject].loc[
            rate_feature_df['state_id']==low_high].loc[
            rate_feature_df['session']==1]
        
        sens_values = subject_df[columns].values[0]
        rate_values_all[low_high]['site1'][i,:] = np.array(np.nan_to_num(sens_values))
    
    for i,row in df_subjects_site2.iterrows():
        subject = row["subject"]
        subject_df = rate_feature_df_site2.loc[
            rate_feature_df_site2['subject']==subject].loc[
            rate_feature_df_site2['state_id']==low_high].loc[
            rate_feature_df_site2['session']==1]
        
        sens_values = subject_df[columns].values[0]
        rate_values_all[low_high]['site2'][i,:] = np.array(np.nan_to_num(sens_values))
    
for low_high in ["low", "high"]:
    p_corrected_list = []
    for sensor in range(0,len(lh_parc)):
        t, p_values = ttest_ind(a=rate_values_all[low_high]['site1'][:,sensor], b=rate_values_all[low_high]['site2'][:,sensor], axis=0)
        _, p_corrected = mne.stats.bonferroni_correction(p_values, alpha=0.05)
        p_corrected_list.append(p_corrected)
    p_values_s1_s2[low_high].append(p_corrected_list)



#################################################################################
###   PLOT SPATIAL DISTRIBUTION OF THE EA VALUES AND STATS BETWEEN THE SITES  ###
#################################################################################

width_mm = 400 # width in mm
height_mm = 60 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

for low_high in ["low", "high"]:
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(width_in, height_in))

    # Mark only differing areas p<0.05 areas
    p_values = 1-np.array(p_values_s1_s2[low_high][0])
    p_values = [1 if x > 0.95 else 0 for x in p_values]
    
    # if there are only zeros in the array, add 0.1 to make the plotting possible
    if sum(p_values) == 0:
        p_values = np.array(p_values)+0.1
        labels_to_stc = mne.labels_to_stc(labels=lh_parc, values=np.nan_to_num(np.array(p_values)), src=src)
    else:
        labels_to_stc = mne.labels_to_stc(labels=lh_parc, values=np.nan_to_num(np.array(p_values)), src=src)

    ss_cont = []

    # Plot the p-values
    for view in ['lat','med']:
        sb_map= sns.color_palette('rocket', as_cmap=True)
        brain = labels_to_stc.plot(subjects_dir=MRI_dir,
                                clim={'kind': 'value' ,'lims':[0.9,0.95,1]},
                                surface='smoothwm',
                                hemi='both',
                                colormap='Blues',
                                background='white',
                                views=view,
                                colorbar=False,
                                transparent=True,
                                smoothing_steps=10)
        
        screenshot = brain.screenshot()
        brain.close()

        nonwhite_pix = (screenshot != 255).any(-1)
        nonwhite_row = nonwhite_pix.any(1)
        nonwhite_col = nonwhite_pix.any(0)
        cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]
        ss_cont.append(cropped_screenshot)

    # Plot also the average rate maps of the wo sites
    plot_limits = []
    for site in ['site1','site2']:
        rates = np.mean(rate_values_all[low_high][site], axis=0)

        labels_to_stc = mne.labels_to_stc(labels=lh_parc, values=np.nan_to_num(np.array(rates)), src=src)

        maxim = np.max(np.nan_to_num(np.array(rates)))
        minim = np.nanmin(np.array(rates))

        clim= {'kind': 'value' ,'lims':[minim,minim+((maxim-(minim))/2),maxim]}
        plot_limits.append(clim)

        for view in ['dorsal']:
            sb_map= sns.color_palette('rocket', as_cmap=True)
            brain = labels_to_stc.plot(subjects_dir=MRI_dir,
                                       clim=clim,
                                       hemi='both',
                                       colormap='magma',
                                       background='white',
                                       views=view,
                                       colorbar=False,
                                       surface='smoothwm')

            screenshot = brain.screenshot()
            brain.close()

            nonwhite_pix = (screenshot != 255).any(-1)
            nonwhite_row = nonwhite_pix.any(1)
            nonwhite_col = nonwhite_pix.any(0)
            cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]
            ss_cont.append(cropped_screenshot)

    # Put the screenshots into the figure and fix axes
    for c in [0,1,2,3]:
        axes[c].imshow(ss_cont[c])
        axes[c].spines['top'].set_visible(False)
        axes[c].spines['right'].set_visible(False)
        axes[c].spines['bottom'].set_visible(False)
        axes[c].spines['left'].set_visible(False)
        axes[c].tick_params(left=False, bottom=False)
        axes[c].set_xticklabels([])
        axes[c].set_yticklabels([])

    for c in [2,3]:
        # Add colorbar below the figure
        ax = axes[c]
        cbar_ax = fig.add_axes([ax.get_position().x0, ax.get_position().y0, ax.get_position().width, 0.02])  
        cbar = mne.viz.plot_brain_colorbar(cbar_ax, plot_limits[c-2], colormap='magma', orientation='horizontal')
        cbar.set_label('FO')

    colors = sns.color_palette("rocket", 5)
    
    # Add box plot
    index_of_precentral = [ 1 if ( 'precentral_12-rh' in lab.name) else 0 for lab in labels_parc].index(1)
    
    df = pd.DataFrame([(group, value[index_of_precentral]) for group, values in rate_values_all[low_high].items() for value in values],columns=['Group', 'Value'])
    sns.boxplot(x='Value', y='Group', data=df, orient='h',ax=axes[4],hue='Group',palette='rocket')
    sns.swarmplot(data=df, x='Value', y='Group', color='grey', size=4,ax=axes[4])

    plt.tight_layout()
    plt.show()
