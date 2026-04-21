"""
Make ICA decomposition, and remove heart artefact, eye blinks and horizontal eye movements.
"""

import mne
import os
import pandas as pd
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
sys.path.append(parent_dir)

from config import (fname, proc)
from settings_hmm_beta import (task, sessions)
from functions_preprocessing import run_ica


# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# For annotations and ICA we use the following frequency limits
lfreq = 1
hfreq = 40

# Iterate over subjects and sessions
for i_ses, session in enumerate(sessions):
    for i_sub, row in df_subjects.iterrows():
        subject = row["subject"]
        
        # Set the raw-file path
        mfilter_path = fname.raw(
            subject=subject,
            ses=session,
            task=task,
            proc=proc,
            )

        if os.path.exists(mfilter_path):  # and not os.path.exists(op_path):

            # Load the raw data
            raw = mne.io.read_raw_fif(mfilter_path).load_data().filter(lfreq,hfreq).resample(sfreq=200)
            
            ica=run_ica('fastica', raw, n_components=50)
            ica.plot_components()
            ica.plot_sources(raw, block=True)

            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=session) + '/ica/'):
                os.makedirs(fname.hmm_bids_dir(subject=subject, ses=session) + '/ica/' )
            
            ica.save( fname.ica(subject=subject,ses=session,task = task ,lfreq = lfreq, hfreq = hfreq), overwrite = True )

