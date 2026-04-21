"""
Create source space and compute forward solution
"""

import os
import pandas as pd
import mne
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from config import (fname, MRI_dir, spacing, N_JOBS, ntri, ico, proc)
from settings_hmm_beta import (task, sessions)

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Iterate over subjects and sessions
for i, row in df_subjects.iterrows():
    subject = row["subject"]
    print('Processing subject: ', subject)

    for ses in sessions:
        print('Processing session: ', ses)

        # Set the raw-file path
        data_path = fname.raw(
            subject=subject,
            ses=ses,
            task=task,
            proc=proc,
        )

        if os.path.exists(data_path):
 
            # Create folder if it does not exist
            if not os.path.exists(fname.hmm_bids_dir(subject=subject, ses=ses) + '/forward/'):  
                os.makedirs(fname.hmm_bids_dir(subject=subject,ses=ses) + '/forward/')
            
            if os.path.exists(fname.src(subject=subject,ses=ses,spacing=spacing)):
                subject_src = mne.read_source_spaces(fname.src(subject=subject, ses=ses,spacing=spacing))
                bem = mne.read_bem_solution(fname.bem(subject=subject, ses=ses,ntri=ntri))
            else:
                # Create source space in individual subject 
                subject_src = mne.setup_source_space(subject=subject,
                                            spacing=spacing,
                                            subjects_dir=MRI_dir,
                                            n_jobs=N_JOBS, add_dist=True)
                mne.write_source_spaces(fname.src(subject=subject,spacing=spacing, ses=ses), subject_src, overwrite=True)

                # Create BEM model
                bem_model = mne.make_bem_model(subject=subject, ico=ico, subjects_dir=MRI_dir,
                                    conductivity=(0.3,))
                if bem_model[0]['ntri'] == ntri:            
                    bem = mne.make_bem_solution(bem_model)
                    mne.write_bem_solution(fname.bem(subject=subject,ntri=bem_model[0]['ntri'], ses=ses),bem,overwrite=True)
                else:
                    raise ValueError('ntri should be %d' % (ntri))

            #Create the forward model, coregistration required (trans_file)
            info = mne.io.read_info(data_path)
            fwd = mne.make_forward_solution(
                info,
                trans=fname.trans(subject= subject, ses=ses, task=task),
                src=subject_src,
                bem=bem,
                meg=True,
                eeg=False,
                mindist=0,
                n_jobs=N_JOBS
                )
            
            # Save forward solution
            mne.write_forward_solution(fname.fwd(subject=subject,ses=ses,task=task, spacing=spacing), fwd,
                            overwrite=True)

