"""
Plot rebound-suppression change of the induced response
"""
import numpy as np
import mne
import seaborn as sns
import pandas as pd
import os
import matplotlib.pyplot as plt

from settings_hmm_beta import (state_mapping, session, group_id)
from config import (MRI_dir, fname, spacing)
from scipy import stats as stats

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])
sorted_map = sorted(state_mapping.items(), key=lambda x: x[1])

# Read the parcels
aparc_sub_parcels = mne.read_labels_from_annot(
    subject = 'fsaverage_sara',
    parc = 'aparc_sub',
    subjects_dir=MRI_dir)

# Read the source space
fname_src = fname.src(hmm_bids_dir=fname.hmm_bids_dir(subject='fsaverage', ses=session), subject='fsaverage', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Read the resting state EA values
EA_feature_df = pd.read_csv(fname.feature_csv(feature = 'EA',job_id=group_id, task='restEO'))


# Let's look at low- an high-beta bands
band_legends = ['low-band','high-band']
rest_labels = ["whole_ts_low", "whole_ts_high"]

induced_great_sum = {band:np.zeros((len(df_subjects), len(aparc_sub_parcels), n_of_tp)) for band in band_legends }
EA_rest_ = {band:np.zeros((len(df_subjects), len(aparc_sub_parcels), 1)) for band in band_legends }

for lh_i, lhfreq in enumerate([(13,20), (20,30)]):

    lfreq = lhfreq[0]
    hfreq = lhfreq[1]

    n_source_points = src[0]['nuse']*2

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

##########################################################
####################    PLOTS   ##########################
##########################################################
# COMPARE ALL 4 STATES TO BL: ONLY SHIFT BASELINE TO THE ZERO

# Set the figure parameters here
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
width_mm = 300 # width in mm
height_mm = 200 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

t_ax = np.linspace(epo_min, epo_max, n_of_tp)
sup_indices = np.where((t_ax >= 0.1) & (t_ax <= 0.4))
reb_indices = np.where((t_ax >= 0.6) & (t_ax <= 1.6))

peak_dif_stcs = []
subject_dif_values = {'low': np.zeros((len(df_subjects), len(lh_parc))),
                      'high': np.zeros((len(df_subjects), len(lh_parc)))}

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(width_in, height_in))

mins_maxs = []
for lh_i, low_high in enumerate(["low-band","high-band"]):
    
    mean_sup = np.mean(induced_great_sum[low_high][:,:,sup_indices][:,:,0,:],axis=2)
    mean_reb = np.mean(induced_great_sum[low_high][:,:,reb_indices][:,:,0,:],axis=2)
    scaled_data = (mean_reb - mean_sup) / EA_rest_[low_high][:,:,0]

    subject_dif_values[low_high] = scaled_data

    plot_data = np.mean(scaled_data,axis=0)


    labels_to_stc = mne.labels_to_stc(labels=lh_parc,
                                      values=np.nan_to_num(np.array(plot_data)),
                                      src=src,)
    
    #contrast_maps.append(labels_to_stc)

    maxim = np.round(np.max(np.nan_to_num(np.array(plot_data))),3)
    minim = np.round(np.nanmin(np.array(plot_data)),3)
    mid = np.round(minim+((maxim-(minim))/2),3)

    mins_maxs.append([minim, maxim])


    clim= {'kind': 'value' ,'lims':[minim,np.round(minim+((maxim-(minim))/2),3),maxim]}

    #clim= {'kind': 'value' ,'lims':[np.percentile(plot_data, 90),
    #                                np.percentile(plot_data, 95),
    #                                np.percentile(plot_data, 100)]}
    
    sb_map= sns.color_palette('rocket', as_cmap=True)
    
    for v_i, view in enumerate(['lateral','medial']):
        #sb_map= sns.color_palette('rocket', as_cmap=True)
        brain = labels_to_stc.plot(subjects_dir=subjects_dir_ave,
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


low_high_band_dif = peak_dif_stcs[0] - peak_dif_stcs[1]

t, p_values = stats.ttest_ind(a=subject_dif_values['low'], b=subject_dif_values['high'], axis=0)
_, p_corrected = mne.stats.bonferroni_correction(p_values, alpha=0.05)

sig_parcel_indices = np.where(p_corrected <= 0.001)[0]
significant_parcels = [lh_parc[parc_i] for parc_i in sig_parcel_indices]

#parcels_lh = np.sum([lh_parcel for lh_parcel in significant_parcels if lh_parcel.hemi == 'lh'])
parcels_rh = np.sum([lh_parcel for lh_parcel in significant_parcels if lh_parcel.hemi == 'rh'])


maxim = np.round(np.max(np.nan_to_num(np.array(mins_maxs))),3)
minim = np.round(np.nanmin(np.array(mins_maxs)),3)
mid = np.round(minim+((maxim-(minim))/2),3)

clim= {'kind': 'value' ,'lims':[0.98*maxim,0.99*maxim,maxim]}

for v_i, view in enumerate(['dorsal']):
    brain = low_high_band_dif.plot(subjects_dir=subjects_dir_ave,
                                hemi='both',
                                clim=clim,
                                views = view,
                                colormap='magma',
                                background='white',
                                colorbar=False,
                                surface='smoothwm')
    
    if parcels_rh != 0:
        brain.add_label(parcels_rh,color=aparc_colors_rh[2],borders=False, alpha=0.8)


    screenshot = brain.screenshot()
    brain.close()

    nonwhite_pix = (screenshot != 255).any(-1)
    nonwhite_row = nonwhite_pix.any(1)
    nonwhite_col = nonwhite_pix.any(0)
    cropped_screenshot = screenshot[nonwhite_row][:, nonwhite_col]

    axes[v_i,2].imshow(cropped_screenshot)


# Add colorbar below the figure
cbar_ax = fig.add_axes([0.2+0.3*2, 0.1, 0.1, 0.02]) # Adjust the position and size as needed
cbar = mne.viz.plot_brain_colorbar(cbar_ax, clim, colormap='magma', orientation='horizontal')


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

plt.savefig('/projects/HMM-beta/HMM_beta_sara/processed/Figures/methods_manu_figures/evoked_gamma_spatial/induced_response_relative_to_rest_task_'+ task  + '.svg')
plt.show()


############################################################

# Set the figure parameters here
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
#plt.rcParams['line.linewidth'] = 0.5
width_mm = 80 # width in mm
height_mm = 55 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Plotting routine
fig, axes = plt.subplots( nrows=1, ncols=1, figsize=(width_in, height_in) )
t_ax = np.linspace(epo_min,epo_max,n_of_tp)
aparc_colors_rh = sns.color_palette("rocket_r", 4).as_hex()
state_idxs = [0,2]
precentral_rh = [1 if par.name == 'precentral_12-rh' else 0 for par in aparc_sub_parcels].index(1)


for state_idx, state in enumerate(['high-band','low-band']):

    curve_data = (induced_great_sum[state][:,parcel_index,:] - EA_rest_[state][:,parcel_index,:] ) / EA_rest_[state][:,parcel_index,:]
    curve = np.mean(curve_data, axis=0)
    print( EA_rest_[state][:,parcel_index,:] )

    axes.plot(t_ax,curve, color=aparc_colors_rh[state_idxs[state_idx]],linewidth=1, label = sorted_state_mapping[state_idxs[state_idx]][0])

axes.set_xlim(epo_min,epo_max)
axes.set_xlabel('Time (s)')
#axes.set_ylabel('Occupancy')
axes.axvline(x=0,linestyle='--', color='black')
axes.axvline(x=0.5,linestyle='--', color='black')
#axes.axhline(y=0,linestyle='--', color='gray')

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('/projects/HMM-beta/HMM_beta_sara/processed/Figures/methods_manu_figures/evoked_gamma_spatial/induced_response_relative_to_rest_task_'+ task  + '_precentral.svg')

plt.show()