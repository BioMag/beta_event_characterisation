"""
Plot rebound-suppression change of the induced response
"""
import numpy as np
import mne
import seaborn as sns
import pandas as pd
import os
import matplotlib.pyplot as plt

from settings_hmm_beta import (group_id, task, task_parameters)
from config import (MRI_dir, fname, spacing)
from scipy import stats as stats

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Set the colors
myColors = sns.color_palette("rocket_r", 4).as_hex()

# Read the parcels
aparc_sub_parcels = mne.read_labels_from_annot(
    subject = 'fsaverage',
    parc = 'aparc_sub',
    subjects_dir=MRI_dir)
lh_parc= [ lab for lab in aparc_sub_parcels if (lab.hemi == 'lh' or lab.hemi == 'rh') ]

# Get the parcels we are interested in
parcel_index = [1 if par.name == 'precentral_12-rh' else 0 for par in aparc_sub_parcels].index(1)

# Read the source space
fname_src = fname.src(megbids_dir=fname.megbids_dir(subject='fsaverage', ses='01'), subject='fsaverage', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Read the resting state EA values
EA_feature_df = pd.read_csv(fname.feature_csv(feature = 'EA',job_id=group_id, task='restEO'))

# We are interested in the low- an high-beta bands
band_legends = ['low-band','high-band']
rest_labels = ["whole_ts_low", "whole_ts_high"]

# Collect the induced responses of the bands from task and rest
induced_great_sum = {band:np.zeros((len(df_subjects), len(aparc_sub_parcels), task_parameters[task]['n_of_tp'])) for band in band_legends }
EA_rest_ = {band:np.zeros((len(df_subjects), len(aparc_sub_parcels), 1)) for band in band_legends }

for lh_i, lhfreq in enumerate([(13,20), (20,30)]):

    lfreq = lhfreq[0]
    hfreq = lhfreq[1]

    for i, row in df_subjects.iterrows():
        subject = row["subject"]
        print(subject)

        stc_path_fname = fname.stc(
            subject = subject,
            ses = '01',
            task = task,
            lfreq = lfreq,
            hfreq=hfreq)
        
        if os.path.exists(stc_path_fname + '-lh.stc'):
            stc_EA = mne.read_source_estimate( stc_path_fname )

            for p_i, parc in enumerate(aparc_sub_parcels):

                label_data = mne.extract_label_time_course(stc_EA,parc, src, mode= 'mean')
                induced_great_sum[band_legends[lh_i]][i,p_i,:] = label_data[0]

                sens_rest_value = EA_feature_df.loc[
                    EA_feature_df['subject']==subject].loc[
                    EA_feature_df['state_id']=='low'].loc[
                    EA_feature_df['session']==1].loc[
                    EA_feature_df['perc']==0].loc[
                    EA_feature_df['measure']==rest_labels[lh_i]]["sens"+str(p_i+1)].values[0]

                EA_rest_[band_legends[lh_i]][i,p_i,0] = sens_rest_value

####################################################################
#################### PLOTS AND STATISTICS ##########################
####################################################################

# Set the figure parameters here
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
width_mm = 300 # width in mm
height_mm = 200 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4


##########################################################
# PLOT THE SPATIAL CHANGE IN THE SUPPRESSION AND REBOUND #
##########################################################


fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(width_in, height_in))
t_ax = np.linspace(task_parameters[task]['epo_tmin'],task_parameters[task]['epo_tmax'],task_parameters[task]['n_of_tp'])

# Collect suppression and rebound mean
sup_indices = np.where((t_ax >= 0.1) & (t_ax <= 0.4))
reb_indices = np.where((t_ax >= 0.6) & (t_ax <= 1.6))

peak_dif_stcs = []
subject_dif_values = {'low': np.zeros((len(df_subjects), len(lh_parc))),
                      'high': np.zeros((len(df_subjects), len(lh_parc)))}

mins_maxs = []

# Loop over the low-beta and high-beta states
for lh_i, low_high in enumerate(["low-band","high-band"]):
    
     # Get the mean over the suppression and rebounds
    mean_sup = np.mean(induced_great_sum[low_high][:,:,sup_indices][:,:,0,:],axis=2)
    mean_reb = np.mean(induced_great_sum[low_high][:,:,reb_indices][:,:,0,:],axis=2)
    scaled_data = (mean_reb - mean_sup) / EA_rest_[low_high][:,:,0]

    subject_dif_values[low_high] = scaled_data

    plot_data = np.mean(scaled_data,axis=0)

    labels_to_stc = mne.labels_to_stc(labels=lh_parc,
                                      values=np.nan_to_num(np.array(plot_data)),
                                      src=src,)

    maxim = np.round(np.max(np.nan_to_num(np.array(plot_data))),3)
    minim = np.round(np.nanmin(np.array(plot_data)),3)
    mid = np.round(minim+((maxim-(minim))/2),3)

    mins_maxs.append([minim, maxim])

    clim= {'kind': 'value' ,'lims':[minim,np.round(minim+((maxim-(minim))/2),3),maxim]}
    
    for v_i, view in enumerate(['lateral','medial']):
        brain = labels_to_stc.plot(subjects_dir=MRI_dir,
                                clim=clim,
                                hemi='both',
                                colormap='magma',
                                background='white',
                                views = view,
                                colorbar=False,
                                surface='smoothwm')
        
        screenshot = brain.screenshot()
        brain.close()

        nonwhite_pix = (screenshot != 255).any(-1)
        nonwhite_row = nonwhite_pix.any(1)
        nonwhite_col = nonwhite_pix.any(0)
        cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]

        axes[v_i,lh_i].imshow(cropped_screenshot)

        # Add colorbar below the figure
        cbar_ax = fig.add_axes([0.2+0.3*lh_i, 0.1, 0.1, 0.02]) # Adjust the position and size as needed
        cbar = mne.viz.plot_brain_colorbar(cbar_ax, clim, colormap='magma', orientation='horizontal')
    
    peak_dif_stcs.append(labels_to_stc)

# statistics and plotting
t, p_values = stats.ttest_ind(a=subject_dif_values['low'], b=subject_dif_values['high'], axis=0)
_, p_corrected = mne.stats.bonferroni_correction(p_values, alpha=0.05)

sig_parcel_indices = np.where(p_corrected <= 0.001)[0]
significant_parcels = [lh_parc[parc_i] for parc_i in sig_parcel_indices]

# Look the right hemisphere
parcels_rh = np.sum([lh_parcel for lh_parcel in significant_parcels if lh_parcel.hemi == 'rh'])

# Plot the significant areas into an empty brain
for v_i, view in enumerate(['dorsal']):
    empty_stc = mne.labels_to_stc(labels=lh_parc, values=np.ones(len(lh_parc)), src=src)

    brain = empty_stc.plot(subjects_dir=MRI_dir,
                            hemi='both',
                            clim={'kind': 'value' ,'lims':[2,3,4]},
                            views = view,
                            colormap='magma',
                            background='white',
                            colorbar=False,
                            surface='smoothwm')
    
    if parcels_rh != 0:
        brain.add_label(parcels_rh,color=myColors[2],borders=False, alpha=0.8)

    screenshot = brain.screenshot()
    brain.close()

    nonwhite_pix = (screenshot != 255).any(-1)
    nonwhite_row = nonwhite_pix.any(1)
    nonwhite_col = nonwhite_pix.any(0)
    cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]

    axes[v_i,2].imshow(cropped_screenshot)

for r in [0,1]:
    for c in [0,1,2]:
        # Hide the axes
        axes[r,c].spines['top'].set_visible(False)
        axes[r,c].spines['right'].set_visible(False)
        axes[r,c].spines['bottom'].set_visible(False)
        axes[r,c].spines['left'].set_visible(False)
        axes[r,c].tick_params(left=False, bottom=False)
        axes[r,c].set_xticklabels([])
        axes[r,c].set_yticklabels([])

plt.show()


##########################################################
# PLOT THE SPATIAL CHANGE IN THE SUPPRESSION AND REBOUND #
##########################################################

# Set the figure parameters here
width_mm = 80 # width in mm
height_mm = 55 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Plotting routine
# Make figure and set the states we are interested in
fig, axes = plt.subplots( nrows=1, ncols=1, figsize=(width_in, height_in) )
t_ax = np.linspace(task_parameters[task]['epo_tmin'],task_parameters[task]['epo_tmax'],task_parameters[task]['n_of_tp'])
band_idxs = [0,2]

for state_idx, band in enumerate(['high-band','low-band']):

    curve_data = (induced_great_sum[band][:,parcel_index,:] - EA_rest_[band][:,parcel_index,:] ) / EA_rest_[band][:,parcel_index,:]
    curve = np.mean(curve_data, axis=0)
    
    # Shift the bl to zero
    curve = curve - np.mean(curve[0:task_parameters[task]['n_of_bl_tp']])
    
    axes.plot(t_ax,curve, color=myColors[band_idxs[state_idx]],linewidth=1, label = band)

axes.set_xlim(task_parameters[task]['epo_tmin'],task_parameters[task]['epo_tmax'])
axes.set_xlabel('Time (s)')
axes.axvline(x=0,linestyle='--', color='black')
axes.axvline(x=0.5,linestyle='--', color='black')

plt.grid()
plt.legend()
plt.tight_layout()
plt.show()