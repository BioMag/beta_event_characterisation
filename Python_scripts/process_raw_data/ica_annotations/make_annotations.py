"""
Annotate data segments that contain muscle artefact.
"""

import mne
import os
import pandas as pd
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
sys.path.append(parent_dir)

from config import (fname, proc)
from settings_hmm_beta import task, sessions


# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# For annotations and ICA we use the following frequency limits
lfreq = 1
hfreq = 40

# Iterate over subjects and sessions
for i, row in df_subjects.iterrows():
    for ses in sessions:
        subject = row["subject"]

        # Set the raw-file path
        mfilter_path = fname.raw(
            subject=subject,
            ses=ses,
            task=task,
            proc=proc,
            )
        
        if os.path.exists(mfilter_path):

            raw = mne.io.read_raw_fif(mfilter_path).load_data().filter(lfreq,hfreq).resample(sfreq=200)

            # Make annotations visually to the raw file
            raw.plot(block=True)

            interactive_annot = raw.annotations

            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=ses) + '/annot/'):
                os.makedirs(fname.hmm_bids_dir(subject=subject, ses=ses) + '/annot/' )

            # Save annotations
            interactive_annot.save( fname.annot(subject=subject,ses=ses,task = task,lfreq = lfreq, hfreq = hfreq), overwrite=True )
