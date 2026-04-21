"""
Make noise covariance.
"""

import pandas as pd
import mne
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from config import fname
from settings_hmm_beta import (task, lfreq, hfreq, sessions)

# Read the subjects from csv
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

for i, row in df_subjects.iterrows():
    subject = row["subject"]
    for ses in sessions:                
        # Get emptyroom data for the noise covariance
        emptyroom_fname = fname.emptyroom(subject=subject,ses=ses)
        
        if os.path.exists(emptyroom_fname):
            raw=mne.io.read_raw_fif(emptyroom_fname).load_data().filter(lfreq,hfreq)
            
            # Calculate noise covariance
            noise_cov = mne.compute_raw_covariance(raw,rank="auto",picks='meg')
            
            # Create folder where the noise_cov is saved, if it does not exist
            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=ses) + '/noise_cov/'):  
                os.makedirs(fname.hmm_bids_dir(subject=subject,ses=ses) + '/noise_cov/')

            # Save noise covariance matrix
            mne.write_cov(fname.noise_cov(subject=subject, ses=ses, task=task, lfreq=lfreq, hfreq=hfreq), noise_cov, overwrite=True)
