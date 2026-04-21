"""
Calculate MNE inverse
"""
import mne
import os
import pandas as pd
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from functions_preprocessing import reject_bad_segments
from config import (fname, spacing, proc, MRI_dir)
from settings_hmm_beta import (task, sessions, lfreq, hfreq, sfreq)
from mne.minimum_norm import make_inverse_operator, apply_inverse_raw
from scipy.io import savemat

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Get labels for FreeSurfer 'aparc_sub' cortical parcellation
labels_parc = mne.read_labels_from_annot(subject='fsaverage', parc='aparc_sub', subjects_dir=MRI_dir)

# Get source space
fname_src = fname.src(subject='fsaverage',ses='01', spacing=spacing)
src = mne.read_source_spaces(fname_src)

# Iterate over subjects and sessions
for i, row in df_subjects.iterrows():
    subject = row["subject"]

    for ses in sessions:
        print('Processing session: ', ses)

        # Set the raw-file path
        data_path = fname.raw(subject=subject,
                                       ses=ses,
                                       task=task,
                                       proc=proc)

        if os.path.exists(data_path):
            
            # Read the forward model
            fwd = mne.read_forward_solution(fname.fwd(subject=subject,ses=ses,task=task, spacing=spacing))
            fwd = mne.convert_forward_solution(fwd, surf_ori=True)

            # Read the noise covariance matrix
            noise_cov = mne.read_cov(fname.noise_cov(subject=subject, ses=ses, task=task, lfreq=lfreq, hfreq=hfreq))
            
            # Read the raw
            raw = mne.io.read_raw_fif(data_path).load_data().filter(1,40)
            
            # Load ica and apply
            ica = mne.preprocessing.read_ica( fname.ica(subject=subject,ses=ses,task = task ,lfreq='1',hfreq='40') )
            ica.apply(raw)

            # Read annotation and reject
            bad_annotation = mne.read_annotations( fname.annot(subject=subject,ses=ses,task = task,lfreq = '1', hfreq = '40') )
            raw.set_annotations(bad_annotation)

            # Then do the rest of the filtering and pick gradiometers
            raw.filter(lfreq,hfreq).resample(sfreq=sfreq).pick_types(meg='grad',stim=True)
            info = raw.info
            
            # Remove the annotations from the data
            raw = reject_bad_segments(raw)

            # Compute inverse operator
            pick_ori = None
            snr = 1.0 
            lambda2 = 1.0 / snr ** 2
            method = "MNE"

            # Create inverse operator
            inverse_operator = make_inverse_operator(info = info, 
                                                forward = fwd, 
                                                noise_cov = noise_cov, 
                                                fixed=True,
                                                loose = 0, 
                                                depth = 0.8,
                                                rank = {'grad':64})

            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=ses) + '/inverse/'):
                os.makedirs(fname.hmm_bids_dir(subject=subject, ses=ses) + '/inverse/' )
            
            # Save inverse operator
            mne.minimum_norm.write_inverse_operator(fname.inv(subject=subject,
                                            ses=ses,
                                            task=task,
                                            spacing = spacing,
                                            lfreq=lfreq,
                                            hfreq=hfreq), inverse_operator, overwrite=True)

            # apply inverse and save source estimate
            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=ses) + '/stcs/'):
                    os.makedirs(fname.hmm_bids_dir(subject=subject, ses=ses) + '/stcs/')

            stc = apply_inverse_raw(raw, inverse_operator, lambda2, method, pick_ori=pick_ori)
            
            stc.save(fname.stc(subject=subject, ses = ses, task=task, lfreq=lfreq, hfreq=hfreq), overwrite=True)

            # Extract the time courses in the label
            label_ts_fsaverage = mne.extract_label_time_course(
                stc, labels_parc, src, mode="pca_flip", allow_empty=True
            )

            # Make the data to modulo 2 => otherwise gives error later
            if label_ts_fsaverage.shape[1] % 2 == 1:
                label_ts_fsaverage = label_ts_fsaverage[:,0:-1]

            # Save data mat
            savemat(fname.data_mat(
                        subject=subject,
                        ses=ses,
                        l_freq=str(lfreq),
                        h_freq=str(hfreq),
                        task=task
                        ),
                        dict(x=label_ts_fsaverage.T),
                    )