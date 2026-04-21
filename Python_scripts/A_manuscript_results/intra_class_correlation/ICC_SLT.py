"""
Calculate and plot intra class correlation of SLT between two sessions 
"""
import pandas as pd
import mne
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from config import (fname, MRI_dir, spacing)
from settings_hmm_beta import (task, group_id)
from pingouin import intraclass_corr

# Read subjects and the characteristics csv
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])
SLT_feature_df = pd.read_csv(fname.feature_csv(feature = 'SLT',job_id=group_id, task=task))

# Read the parcels
labels_parc = mne.read_labels_from_annot(subject = 'fsaverage', parc = 'aparc_sub', subjects_dir=MRI_dir)
lh_parc= [ lab for lab in labels_parc if (lab.hemi == 'lh' or lab.hemi == 'rh') ]

# Find the label index and object of the pre-selected pre-central parcel
peak_label = ['precentral_12-rh',]
peak_label_idx = []
peak_label_object = []

for nl_i, nl in enumerate(peak_label):
    index = [ 1 if ( nl in lab.name) else 0 for lab in labels_parc].index(1)
    peak_label_idx.append(index)

    lh_labels = [ lab for lab in labels_parc if ( nl in lab.name) ]
    peak_label_object.append(lh_labels[0])

# Read the source space
fname_src = fname.src(hmm_bids_dir=fname.hmm_bids_dir(subject='fsaverage', ses='01'), subject='fsaverage', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Collect the ICC values from each parcel and SLT values from the pre-central parcel
parcel_values = {'low':[],'high':[]}
peak_sens_vals = {'low':{'ses01':[], 'ses02':[]},'high':{'ses01':[], 'ses02':[]}}

# Loop over all parcels to get ICC from each
for sensor in range(0,len(lh_parc)):
    sens_id = 'sens' + str(sensor+1)
    
    for low_high in ["low", "high"]:
        # Make empty list for SLT values
        char_values = {'subjects':[],'session':[],'score':[]}
        
        for ses in [1,2]:
            for i,row in df_subjects.iterrows():
                subject = row["subject"]
                subject_df = SLT_feature_df.loc[
                    SLT_feature_df['subject']==subject].loc[
                    SLT_feature_df['state_id']==low_high].loc[
                    SLT_feature_df['session']==ses].loc[
                    SLT_feature_df['measure']=='mean']
                sens_values = subject_df[sens_id].values[0]*5
                
                char_values['subjects'].append(subject)
                char_values['session'].append(ses)
                char_values['score'].append(sens_values)
                
                if sensor == peak_label_idx[0]:
                    peak_sens_vals[low_high]['ses0'+str(ses)].append(sens_values)

        # Calculate the ICC
        df = pd.DataFrame(data=char_values)
        ICC_parcel = intraclass_corr(data=df, targets='subjects', raters='session', ratings='score',nan_policy='omit').round(3)
        parcel_values[low_high].append(ICC_parcel['ICC'][2])

###############################################
###   PLOT SPATIAL DISTRIBUTION OF THE ICC  ###
###############################################

# Set the figure parameters here
width_mm = 250 # width in mm
height_mm = 60 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Plot for low-beta state and high-beta state
for low_high in ["low", "high"]:

    # Create figure
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(width_in, height_in))
    cmap = sns.color_palette("Blues", as_cmap=True)
    
    # Find the ICC values of the state
    ICC_values = parcel_values[low_high]

    # Mark different color for different ICC limits
    ICC_values = [10 if i>=0.75 else i for i in ICC_values]
    ICC_values = [7 if 0.6<=i<0.75 else i for i in ICC_values]
    ICC_values = [4 if 0.4<=i<0.6 else i for i in ICC_values]
    ICC_values = [1 if i<0.4 else i for i in ICC_values]
    
    # Make sourceEstimate object of the ICC values
    labels_to_stc = mne.labels_to_stc(labels=lh_parc, values=np.nan_to_num(np.array(ICC_values)), src=src)

    # Collect screenshots of different views
    ss_cont = []
    for view in ['lat','med','dorsal']:
        sb_map= sns.color_palette('rocket', as_cmap=True)
        brain = labels_to_stc.plot(subjects_dir=MRI_dir,
                                clim={'kind': 'value' ,'lims':[1,5,10]},
                                surface='smoothwm',
                                hemi='both',
                                colormap=cmap,
                                background='white',
                                views=view,
                                colorbar=False,
                                transparent=False,
                                smoothing_steps=5)
        
        # Put borders of the parcels
        [brain.add_label(lab,color='lightgrey', alpha=1, borders=0.1) for lab in labels_parc]

        screenshot = brain.screenshot()
        brain.close()

        nonwhite_pix = (screenshot != 255).any(-1)
        nonwhite_row = nonwhite_pix.any(1)
        nonwhite_col = nonwhite_pix.any(0)
        cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]
        ss_cont.append(cropped_screenshot)
    
    for c in [0,1,2]:
        axes[c].imshow(ss_cont[c])
        axes[c].spines['top'].set_visible(False)
        axes[c].spines['right'].set_visible(False)
        axes[c].spines['bottom'].set_visible(False)
        axes[c].spines['left'].set_visible(False)
        axes[c].tick_params(left=False, bottom=False)
        axes[c].set_xticklabels([])
        axes[c].set_yticklabels([])
    
plt.show()

##########################################################
###   PLOT THE FO VALUES IN THE RH PRE-CENTRAL PARCEL  ###
##########################################################

# Set the figure parameters here
width_mm = 60 # width in mm
height_mm = 50 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

x_lims = [(100,225),(120,260)]

# Plot separately for low- and high-beta states
for lh_i, low_high in enumerate(["low", "high"]):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(width_in, height_in))

    colors = sns.color_palette("rocket", 1)
    ICC_values = parcel_values[low_high]
    ICC_top = ICC_values[peak_label_idx[0]]

    # Make scatter plot
    axes.scatter(peak_sens_vals[low_high]['ses01'],
                peak_sens_vals[low_high]['ses02'],
                color=colors[0],
                marker='o',
                label = 'ICC - ' + str(ICC_top))

    # Set axes labels and limits
    axes.set_xlabel('session 1')
    axes.set_ylabel('session 2')
    axes.set_xlim(x_lims[lh_i][0], x_lims[lh_i][1])
    axes.set_ylim(x_lims[lh_i][0], x_lims[lh_i][1])
    axes.legend(loc='lower right')
    plt.tight_layout()

plt.show()
