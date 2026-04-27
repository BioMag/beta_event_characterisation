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
study_path = '/projects/HMM-beta/data/'

# Source spacing
spacing = "ico4"
ntri=5120
ico=4

N_JOBS = 4

proc = 'tsss'


###############################################################################
# Folders (TODO: decide on folder structure)

# datapath
data_path = os.path.join(study_path, "MEG_derivatives")

# path for the processed data (in HMM-MAR)
processed_dir = os.path.join(study_path, "HMM_derivatives/")

# Scaled MRÍ director
MRI_dir = os.path.join(study_path, 'MRI/')

###############################################################################
# Templates for filenames
fname = FileNames()

# Some directories<<<<<<<<<<<<<<<<<<<<<<<<<
fname.add("study_path", study_path)
fname.add("data_path", data_path)
fname.add("processed_dir", processed_dir)

fname.add("megbids_dir", "{data_path}/{subject}/ses-{ses}/meg")
fname.add("hmm_dir", "{processed_dir}/{subject}/")


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
    "ica", "{megbids_dir}/ica/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-ica.fif"
)
# Annotations
fname.add(
    'annot', "{megbids_dir}/annot/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-annot.fif"
)
# Source level files
# Noise covariance
fname.add(
    "noise_cov", "{megbids_dir}/noise_cov/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-cov.fif"
)
# source space
fname.add(
    'src', '{megbids_dir}/forward/{subject}-{spacing}-src.fif'
)
# BEm solution
fname.add(
    'bem', '{megbids_dir}/forward/{subject}-{ntri}-bem-sol.fif'
)
#Forward model
fname.add(
    'fwd', '{megbids_dir}/forward/{subject}_ses-{ses}_task-{task}_{spacing}-fwd.fif'
)
# inverse solution
fname.add(
    'inv', '{megbids_dir}/inverse/{subject}_ses-{ses}_task-{task}-{spacing}_lfreq-{lfreq}-hfreq-{hfreq}-inv.fif'
)
# Source estimate
fname.add(
    'stc', '{megbids_dir}/stcs/{subject}_ses-{ses}_task-{task}_lfreq-{lfreq}-hfreq-{hfreq}-stc'
)
# Save .mat file
fname.add(
    "data_mat",
    "{data_path}/{subject}/ses-{ses}/data_mat/{subject}_lfreq-{l_freq}_hfreq-{h_freq}_task-{task}_mat.mat",
)


# Source estimate characteristics
# hmm-dual
fname.add(
    "hmm_dual", "{hmm_dir}/dual_hmm_ses_{ses}_{job_id}_task_{task}.mat"
)
fname.add(
    "FO_dual", "{hmm_dir}/FO_ses_{ses}_{job_id}_task_{task}.mat"
)
fname.add(
    "SLT_dual", "{hmm_dir}/SLT_parameters_ses_{ses}_{job_id}_task_{task}.mat"
)
fname.add(
    "TPM_dual", "{hmm_dir}/TPM_parameters_ses_{ses}_{job_id}_task_{task}.mat"
)

# Feature CSVs
fname.add(
    "feature_csv", "{study_path}/characteristics_csvs/{feature}_csv/feature-{feature}_task-{task}_jobid-{job_id}.csv"
)

