"""
===========
Config file
===========

Configuration parameters for the study.
"""

import os

from fnames import FileNames

user = os.environ["USER"]


# path to where the data exists
#study_path = '/data/path/'
study_path = '/projects/HMM-beta/HMM_beta_sara/test_data'

# Source spacing
spacing = "ico4"
ntri=5120
ico=4

N_JOBS = 4

proc = 'tsss'


###############################################################################
# Folders (TODO: decide on folder structure)

# scimeg datapath
data_path = os.path.join(study_path, "BIDS")

# path for the processed data (in HMM-MAR)
processed_dir = os.path.join(study_path, "processed/")

# Scaled MRÍ director
MRI_dir = os.path.join(study_path, 'MRI/')

###############################################################################
# Templates for filenames
fname = FileNames()

# Some directories<<<<<<<<<<<<<<<<<<<<<<<<<
fname.add("data_path", data_path)
fname.add("processed_dir", processed_dir)

fname.add("megbids_dir", "{data_path}/{subject}/ses-{ses}/meg")
fname.add("hmm_bids_dir", "{processed_dir}/{subject}/ses-{ses}/")


# Sensor-level files
# Raw files in HMM folders
fname.add(
    "raw", "{megbids_dir}/{subject}_ses-{ses}_task-{task}_proc-{proc}.fif"
)
# Emptyroom file for the noise covariance
fname.add(
    "emptyroom", "{megbids_dir}/emptyroom_tsss.fif"
)
# Trans-file (coregistration)
fname.add('trans',
    '{megbids_dir}/{subject}_ses-{ses}_{task}-coreg-trans.fif'
)



# ica decompositions
fname.add(
    "ica", "{hmm_bids_dir}/ica/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-ica.fif"
)
# Annotations
fname.add(
    'annot', "{hmm_bids_dir}/annot/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-annot.fif"
)

# Source level files
# Noise covariance
fname.add(
    "noise_cov", "{hmm_bids_dir}/noise_cov/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-cov.fif"
)
# source space
fname.add(
    'src', '{hmm_bids_dir}/forward/{subject}-{spacing}-src.fif'
)
# BEm solution
fname.add(
    'bem', '{hmm_bids_dir}/forward/{subject}-{ntri}-bem-sol.fif'
)
#Forward model
fname.add(
    'fwd', '{hmm_bids_dir}/forward/{subject}_ses-{ses}_task-{task}_{spacing}-fwd.fif'
)
# inverse solution
fname.add(
    'inv', '{hmm_bids_dir}/inverse/{subject}_ses-{ses}_task-{task}-{spacing}_lfreq-{lfreq}-hfreq-{hfreq}-inv.fif'
)
# Source estimate
fname.add(
    'stc', '{hmm_bids_dir}/stcs/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-stc'
)
# Save .mat file
fname.add(
    "data_mat",
    "{hmm_bids_dir}/data_mat/{subject}_lfreq-{l_freq}_hfreq-{h_freq}_task-{task}_mat.mat",
)


# Source estimate characteristics
fname.add(
    'stc_char', '{hmm_bids_dir}/stcs_char/{subject}_ses-{ses}_task-{task}_char-{char}_stc'
)

# hmm-dual
fname.add(
    "hmm_dual", "{processed_dir}/{subject}/ses-{ses}/HMM_output/dual_hmm_{job_id}_task_{task}.mat"
)
fname.add(
    "FO_dual", "{processed_dir}/{subject}/ses-{ses}/HMM_output/FO_{job_id}_task_{task}.mat"
)
fname.add(
    "SLT_dual", "{processed_dir}/{subject}/ses-{ses}/HMM_output/SLT_parameters_{job_id}_task_{task}.mat"
)
fname.add(
    "TPM_dual", "{processed_dir}/{subject}/ses-{ses}/HMM_output/TPM_parameters_{job_id}_task_{task}.mat"
)


# Feature CSVs
fname.add(
    "feature_csv", "{processed_dir}/characteristics_csvs/{feature}_csv/feature-{feature}_task-{task}_jobid-{job_id}.csv"
)

# statistical
fname.add(
    "coincidence_areas", "{processed_dir}/coincidences/coincidence_maps/statistical-map_state-{state}_side-effect-{side_effect}_limitperc-{perc}.stc"
)
fname.add(
    "stat_clust", "{processed_dir}/spatial_clusters/clusters_state-{state}_task-{task}"
)
fname.add(
    "coincidence", "{processed_dir}/coincidences/all_to_all_tp_similarity_task-{task}.pkl"
)
