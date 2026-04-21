import scipy
import numpy as np
import matplotlib.pyplot as plt
import mne
import seaborn as sns
import pandas as pd

from config import (subjects_dir_ave, fname, spacing)
from matplotlib import font_manager
from scipy import stats as stats
from settings_hmm_beta import state_mapping



# Add font
font_dirs = ["/usr/share/fonts/truetype/msttcorefonts/"]  # The path to the custom font file.
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    font_manager.fontManager.addfont(font_file)

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/scimeg_controls_all.txt", names=["subject"])

# Set the colors and task
myColors = sns.color_palette("rocket_r", 4).as_hex()
task = 'leftCKCevoked' #'leftmotor' #'leftCKCevoked'
n_of_tp = 701 #701 #1401 #701
n_of_bl_tp = 100 #88 #192
epo_min = -0.5 #-0.5
epo_max = 3 #3
job_id_sci = "_job_source_2"

sorted_state_mapping = sorted(state_mapping.items(), key=lambda item: item[1])

# Read the parcels
aparc_sub_parcels = mne.read_labels_from_annot(
    subject = 'fsaverage_sara',
    parc = 'aparc_sub',
    subjects_dir=subjects_dir_ave)

lh_parc= [ lab for lab in aparc_sub_parcels if (lab.hemi == 'lh' or lab.hemi == 'rh') ]

parcel_index = [1 if par.name == 'precentral_12-rh' else 0 for par in aparc_sub_parcels].index(1)
    
fname_src = fname.src(hmm_bids_dir=fname.hmm_bids_dir(subject='fsaverage_sara', ses='01'), subject='fsaverage_sara', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Collect the gamma-paths from the states
gamp_great_sum = {state[0]:np.zeros((len(df_subjects), len(aparc_sub_parcels), n_of_tp)) for state in sorted_state_mapping }

# Read the resting state FO values
FO_feature_df = pd.read_csv(fname.feature_csv(group_or_single = 'group', feature = 'FO',job_id=job_id_sci, task='restEO'))

FO_rest_ = {state[0]:np.zeros((len(df_subjects), len(aparc_sub_parcels), 1)) for state in sorted_state_mapping }


for i, row in df_subjects.iterrows():
    subject = row["subject"]
    #fig, axes = plt.subplots(nrows=len(labels), ncols=1, figsize=(5, 6))
    print(subject)
    
    gamma_path_fname = fname.hmm_dual(
        subject = subject,
        group_or_single="group",
        sensor_type="sensors_concat",
        session = '01',
        job_id = job_id_sci,
        task = task)
    
    gamma_path_rest_fname = fname.hmm_dual(
        subject = subject,
        group_or_single="group",
        sensor_type="sensors_concat",
        session = '01',
        job_id = job_id_sci,
        task = 'restEO')

    hmm_dual_dict = scipy.io.loadmat( gamma_path_fname )
    hmm_dual_rest_dict = scipy.io.loadmat( gamma_path_rest_fname )
    
    for i_s, s_id in enumerate(aparc_sub_parcels):

        gam_path = hmm_dual_dict["output"]["sens"+str(i_s+1)][0][0]["gamma"][0][0]

        gam_path_rest = hmm_dual_rest_dict["output"]["sens"+str(i_s+1)][0][0]["gamma"][0][0]

        # Add 8 zeros to the front and after so that all of the epochs are the same length
        gam_path = gam_path[(n_of_tp-8):-(n_of_tp-8),:]

        n_epochs = int(gam_path.shape[0]/n_of_tp)

        f_i = 0
        l_i = n_of_tp

        sum_gampath = np.zeros((n_of_tp,4))

        for _ in range(0,n_epochs):
            sum_gampath += gam_path[f_i:l_i,:]
            f_i += n_of_tp
            l_i += n_of_tp

        for state_i, state in enumerate(sorted_state_mapping):
            gamp_great_sum[state[0]][i,i_s,:] = sum_gampath[:,state_i] / n_epochs

            #sens_rest_value = FO_feature_df.loc[
            #    FO_feature_df['subject']==subject].loc[
            #    FO_feature_df['state_id']==state[0]].loc[
            #    FO_feature_df['session']==1]["sens"+str(i_s+1)].values[0]
            
            sens_rest_value = np.mean(gam_path_rest[:,state_i])

            FO_rest_[state[0]][i,i_s,0] = sens_rest_value


##########################################################
#################### STATISTICS ##########################
##########################################################
# COMPARE ALL 4 STATES TO BL: ONLY SHIFT BASELINE TO THE ZERO

# Set the figure parameters here
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
width_mm = 200 # width in mm
height_mm = 150 # height in mm
width_in = width_mm / 25.4
height_in = height_mm / 25.4

t_ax = np.linspace(epo_min, epo_max, n_of_tp)
sup_indices = np.where((t_ax >= 0.1) & (t_ax <= 0.4))
reb_indices = np.where((t_ax >= 0.6) & (t_ax <= 1.6))

peak_dif_stcs = []
subject_dif_values = {'low': np.zeros((len(df_subjects), len(lh_parc))),
                      'high': np.zeros((len(df_subjects), len(lh_parc)))}

# Then create one figure for all plots (finally add boxplot from certain label)
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(width_in, height_in))

mins_maxs = []

for lh_i, low_high in enumerate(["low","high"]):
    
    mean_sup = np.mean(gamp_great_sum[low_high][:,:,sup_indices][:,:,0,:],axis=2)
    mean_reb = np.mean(gamp_great_sum[low_high][:,:,reb_indices][:,:,0,:],axis=2)
    scaled_data = (mean_reb - mean_sup) / FO_rest_[low_high][:,:,0]

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
    
    #sb_map= sns.color_palette('rocket', as_cmap=True)
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

low_high_state_dif = peak_dif_stcs[0] - peak_dif_stcs[1]

t, p_values = stats.ttest_ind(a=subject_dif_values['low'], b=subject_dif_values['high'], axis=0)
_, p_corrected = mne.stats.bonferroni_correction(p_values, alpha=0.05)

sig_parcel_indices = np.where(p_corrected <= 0.001)[0]
significant_parcels = [lh_parc[parc_i] for parc_i in sig_parcel_indices]

#parcels_lh = np.sum([lh_parcel for lh_parcel in significant_parcels if lh_parcel.hemi == 'lh'])
parcels_rh = np.sum([lh_parcel for lh_parcel in significant_parcels if lh_parcel.hemi == 'rh'])

#significant_stc = low_high.in_label(parcels_rh)

#clim= {'kind': 'value' ,'lims':[np.percentile(low_high.data, 90),
#                                np.percentile(low_high.data, 95),
#                                np.percentile(low_high.data, 100)]}

maxim = np.round(np.max(np.nan_to_num(np.array(mins_maxs))),3)
minim = np.round(np.nanmin(np.array(mins_maxs)),3)
mid = np.round(minim+((maxim-(minim))/2),3)

clim= {'kind': 'value' ,'lims':[0.98*maxim,0.99*maxim,maxim]}
aparc_colors_rh = sns.color_palette("rocket_r", 4).as_hex()

for v_i, view in enumerate(['dorsal']):
    brain = low_high_state_dif.plot(subjects_dir=subjects_dir_ave,
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

plt.savefig('/projects/HMM-beta/HMM_beta_sara/processed/Figures/methods_manu_figures/evoked_gamma_spatial/gamma_relative_to_rest_task_'+ task  + '.svg')
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

for state_idx, state in enumerate(['high','low']):

    curve_data = (gamp_great_sum[state][:,parcel_index,:] - FO_rest_[state][:,parcel_index,:] ) / FO_rest_[state][:,parcel_index,:]
    curve = np.mean(curve_data, axis=0)
    
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
plt.savefig('/projects/HMM-beta/HMM_beta_sara/processed/Figures/methods_manu_figures/evoked_gamma_spatial/gamma_relative_to_rest_task_'+ task  + '_precentral.svg')

plt.show()