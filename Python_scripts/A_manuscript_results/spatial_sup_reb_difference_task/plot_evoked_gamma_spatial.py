"""
Plot rebound-suppression change of the evoked state probabilities
"""
import scipy
import numpy as np
import matplotlib.pyplot as plt
import mne
import seaborn as sns
import pandas as pd

from config import (MRI_dir, fname, spacing)
from scipy import stats as stats
from settings_hmm_beta import (state_mapping, task_parameters, task, group_id, session)


# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Set the colors
myColors = sns.color_palette("rocket_r", 4).as_hex()

# Sor the states
sorted_state_mapping = sorted(state_mapping.items(), key=lambda item: item[1])

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

# Collect the gamma-paths of the states from task and rest
gamp_great_sum = {state[0]:np.zeros((len(df_subjects), len(aparc_sub_parcels), task_parameters[task]['n_of_tp'])) for state in sorted_state_mapping }
gam_path_mean_rest_ = {state[0]:np.zeros((len(df_subjects), len(aparc_sub_parcels), 1)) for state in sorted_state_mapping }

# Loop over the subject 
for i, row in df_subjects.iterrows():
    subject = row["subject"]
    print(subject)
    
    # Set gamma-path names for task and rest
    gamma_path_fname = fname.hmm_dual(
        subject = subject,
        ses = session,
        job_id = group_id,
        task = task)
    
    gamma_path_rest_fname = fname.hmm_dual(
        subject = subject,
        ses = session,
        job_id = group_id,
        task = 'restEO')

    # Load the dual-HMMs to dictionaries
    hmm_dual_dict = scipy.io.loadmat( gamma_path_fname )
    hmm_dual_rest_dict = scipy.io.loadmat( gamma_path_rest_fname )
    
    # Loop over the parcels
    for i_s, s_id in enumerate(aparc_sub_parcels):
        
        gam_path = hmm_dual_dict["output"]["sens"+str(i_s+1)][0][0]["gamma"][0][0]

        gam_path_rest = hmm_dual_rest_dict["output"]["sens"+str(i_s+1)][0][0]["gamma"][0][0]

        # Take the first and last epochs away
        gam_path = gam_path[(task_parameters[task]['n_of_tp']-8):-(task_parameters[task]['n_of_tp']-8),:]

        n_epochs = int(gam_path.shape[0]/task_parameters[task]['n_of_tp'])

        f_i = 0
        l_i = task_parameters[task]['n_of_tp']

        sum_gampath = np.zeros((task_parameters[task]['n_of_tp'],4))

        for _ in range(0,n_epochs):
            sum_gampath += gam_path[f_i:l_i,:]
            f_i += task_parameters[task]['n_of_tp']
            l_i += task_parameters[task]['n_of_tp']

        for state_i, state in enumerate(sorted_state_mapping):
            gamp_great_sum[state[0]][i,i_s,:] = sum_gampath[:,state_i] / n_epochs

            sens_rest_value = np.mean(gam_path_rest[:,state_i])
            gam_path_mean_rest_[state[0]][i,i_s,0] = sens_rest_value


####################################################################
#################### PLOTS AND STATISTICS ##########################
####################################################################

# Set the figure parameters here
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
width_mm = 200 # width in mm
height_mm = 150 # height in mm
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
for lh_i, low_high in enumerate(["low","high"]):
    
    # Get the mean over the suppression and rebounds
    mean_sup = np.mean(gamp_great_sum[low_high][:,:,sup_indices][:,:,0,:],axis=2)
    mean_reb = np.mean(gamp_great_sum[low_high][:,:,reb_indices][:,:,0,:],axis=2)
    scaled_data = (mean_reb - mean_sup) / gam_path_mean_rest_[low_high][:,:,0]

    subject_dif_values[low_high] = scaled_data

    plot_data = np.mean(scaled_data,axis=0)

    labels_to_stc = mne.labels_to_stc(labels=lh_parc,
                                      values=np.nan_to_num(np.array(plot_data)),
                                      src=src,)
    
    
    maxim = np.round(np.max(np.nan_to_num(np.array(plot_data))),3)
    minim = np.round(np.nanmin(np.array(plot_data)),3)
    mid = np.round(minim+((maxim-(minim))/2),3)

    mins_maxs.append([minim, maxim])

    clim = {'kind': 'value' ,'lims':[minim,np.round(minim+((maxim-(minim))/2),3),maxim]}

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

# Make figure and set the states we are interested in
fig, axes = plt.subplots( nrows=1, ncols=1, figsize=(width_in, height_in) )
t_ax = np.linspace(task_parameters[task]['epo_tmin'],task_parameters[task]['epo_tmax'],task_parameters[task]['n_of_tp'])
state_idxs = [0,2]

for state_idx, state in enumerate(['high','low']):

    curve_data = (gamp_great_sum[state][:,parcel_index,:] - gam_path_mean_rest_[state][:,parcel_index,:] ) / gam_path_mean_rest_[state][:,parcel_index,:]
    curve = np.mean(curve_data, axis=0)

    # Shift the bl to zero
    curve = curve - np.mean(curve[0:task_parameters[task]['n_of_bl_tp']])
    
    axes.plot(t_ax,curve, color=myColors[state_idxs[state_idx]],linewidth=1, label = sorted_state_mapping[state_idxs[state_idx]][0])

axes.set_xlim(task_parameters[task]['epo_tmin'],task_parameters[task]['epo_tmax'])
axes.set_xlabel('Time (s)')
axes.axvline(x=0,linestyle='--', color='black')
axes.axvline(x=0.5,linestyle='--', color='black')

plt.grid()
plt.legend()
plt.tight_layout()
plt.show()










