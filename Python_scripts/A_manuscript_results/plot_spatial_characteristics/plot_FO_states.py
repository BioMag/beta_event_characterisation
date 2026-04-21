"""
Plot FO values across the whole brain
"""
import pandas as pd
import mne
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from config import (fname, MRI_dir,spacing)
from settings_hmm_beta import (task, group_id, sessions)

# Read subjects and the characteristics csv
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])
FO_feature_df = pd.read_csv(fname.feature_csv(feature = 'FO',job_id=group_id, task=task))
session = sessions[0]

# Read the parcels
parcels = mne.read_labels_from_annot(subject = 'fsaverage', parc = 'aparc_sub', subjects_dir=MRI_dir)

# Read the source space
fname_src = fname.src(hmm_bids_dir=fname.hmm_bids_dir(subject='fsaverage', ses=session), subject='fsaverage', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Make columns to take from the df
columns = []
for sens in range(0,len(parcels)):
    columns.append('sens' + str(sens+1))

# Loop over the states you want to plot
for low_high in ["low", "high", "gam", "bl"]:
    FO_values_all = np.zeros(len(parcels))

    for i,row in df_subjects.iterrows():
        subject = row["subject"]
        subject_df = FO_feature_df.loc[
            FO_feature_df['subject']==subject].loc[
            FO_feature_df['state_id']==low_high].loc[
            FO_feature_df['session']==int(session[1])]

        sens_values = subject_df[columns].values[0]
        FO_values_all += np.array(sens_values)

    FO_values_all = FO_values_all/len(df_subjects)
    
    # Put the parcel values to source estimate
    labels_to_stc = mne.labels_to_stc(labels=parcels, values=np.nan_to_num(np.array(FO_values_all)), src=src)

    # Create limits for the plot
    maxim = np.max(np.nan_to_num(np.array(FO_values_all)))
    minim = np.nanmin(np.array(FO_values_all))
    clim= {'kind': 'value' ,'lims':[minim,minim+((maxim-(minim))/2),maxim]}

    ss_cont = []
    for view in ['lateral','medial','dorsal']:
        sb_map= sns.color_palette('rocket', as_cmap=True)
        brain = labels_to_stc.plot(subjects_dir=MRI_dir,
            clim=clim, hemi='both',
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

    # Then create one figure for all plots
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 6))
            
    im1 = axes[0].imshow(ss_cont[0])
    axes[1].imshow(ss_cont[1])
    axes[2].imshow(ss_cont[2])

    for r in [0,1,2]:
        #for c in [0,1]:
        # Hide the axes
        axes[r].spines['top'].set_visible(False)
        axes[r].spines['right'].set_visible(False)
        axes[r].spines['bottom'].set_visible(False)
        axes[r].spines['left'].set_visible(False)
        axes[r].tick_params(left=False, bottom=False)
        axes[r].set_xticklabels([])
        axes[r].set_yticklabels([])

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.7, 0.25, 0.02, 0.5])
    cbar = mne.viz.plot_brain_colorbar(cbar_ax, clim,colormap='magma')
    cbar.set_label('FO')

plt.show()
