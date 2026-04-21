"""
Make CSV file of the event amplitudes.
CSV has (subjects as rows and parcels values in the columns)
"""
import numpy as np
import scipy
import pandas as pd
import os
import scipy
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from config import fname
from functions_make_characteristics import beta_amplitude_envelope, amplitude_envelope, envelope_cutted_list
from settings_hmm_beta import (state_mapping, lfreq, hfreq,
                               n_labels, lag, task, group_id, sessions, sfreq)

df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Define the amplitude measures that are wanted to be extracted
measures = ["mean", "median", "std", "max", "robust_max"]
cut_percentage = [0, 75] # Set percentage limit for the amplitude

# Sorted mapping
sorted_map = sorted(state_mapping.items(), key=lambda x: x[1])

# Make columns for the csv
columns = ["subject","session","state_id","measure"]
for sens in range(0,n_labels):
    columns.append('sens' + str(sens+1))

# make empty dataframe
save_dict = pd.DataFrame(columns=columns)

# Loop over sessions
for j, session in enumerate(sessions):
    # Loop over subjects
    for i, row in df_subjects.iterrows():
        subject = row["subject"]
        print(subject)

        hmm_dual_path = fname.hmm_dual(
            subject = subject,
            ses = session,
            job_id = group_id,
            task = task)
        
        data_path = fname.data_mat(
            subject=subject,
            ses=session,
            l_freq=str(lfreq),
            h_freq=str(hfreq),
            task = task)

        # CHeck if both data file and hmm dual exists
        if os.path.exists(hmm_dual_path) and os.path.exists(data_path):
            print("found dual")
            # Load the data
            data = scipy.io.loadmat(data_path)
            data = data["x"].T
            data = data[:,lag:-lag]

            # Load the dual
            hmm_dual_dict = scipy.io.loadmat( hmm_dual_path )
            length = len(hmm_dual_dict["output"]["sens1"][0][0]["viterbi"][0][0])

            # Time frequency transform of data.
            power, freqs = beta_amplitude_envelope(
                data, sfreq, lower_freq=lfreq, upper_freq=hfreq
            )

            # Make amplitude envelope for the specified frequency range
            # This results into (n_channels, n_timepoints) shaped array
            # Meaning one amplitude envelope for one channel
            amp_env = amplitude_envelope(lfreq, hfreq, freqs, power)
            amp_env_low = amplitude_envelope(13, 20, freqs, power)
            amp_env_high = amplitude_envelope(20, 30, freqs, power)

            # Calculate the amplitude characteristics for all states separately
            for s_key_i, s_key in enumerate(sorted_map):
                # Loop over the different percentages where to cut the AE
                for p, perc in enumerate(cut_percentage):
                    dict_append_mean = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'mean','perc':perc}
                    dict_append_median = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'median','perc':perc}
                    dict_append_std = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'std','perc':perc}
                    dict_append_max = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'max','perc':perc}
                    dict_append_robust_max = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'robust_max','perc':perc}
                    dict_append_whole_ts = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'whole_ts','perc':perc}
                    dict_append_whole_ts_low = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'whole_ts_low','perc':perc}
                    dict_append_whole_ts_high = {'subject': subject,'session':session,'state_id':s_key[0],'measure':'whole_ts_high','perc':perc}
     
                    for sensor in range(0, n_labels):
                        # Make the sensor id and take the viterbi-path
                        sens_id = "sens" + str(sensor + 1)
                        vmap_sens = hmm_dual_dict["output"][sens_id][0][0]["viterbi"][0][0]
                        
                        # Take the amplitude envelope of the label
                        ae_sens = amp_env[sensor, :]
                        ae_sens_low = amp_env[sensor, :]
                        ae_sens_high = amp_env[sensor, :]

                        # If the data and amplitude envelope are not equal length,
                        # take corresponding part of the vmap
                        if len(amp_env[0, :]) != len(data[0, :]):
                            len_dif = len(data[0, :]) - len(amp_env[0, :])
                            vmap_sens = vmap_sens[0 : -len_dif]
                        
                        #Calculate whole time series mean to the dictionary
                        mean_of_whole_ts = np.mean( ae_sens[ ae_sens >=  np.quantile(ae_sens, perc/100)] )
                        dict_append_whole_ts[sens_id] = mean_of_whole_ts

                        mean_of_whole_ts_low = np.mean( ae_sens_low[ ae_sens_low >=  np.quantile(ae_sens_low, perc/100)] )
                        dict_append_whole_ts_low[sens_id] = mean_of_whole_ts_low

                        mean_of_whole_ts_high = np.mean( ae_sens_high[ ae_sens_high >=  np.quantile(ae_sens_high, perc/100)] )
                        dict_append_whole_ts_high[sens_id] = mean_of_whole_ts_high
                        
                        # Make "zeroed envelope", which means that all other values
                        # are zeros but the ones that are in the state
                        vmap_multiply = (vmap_sens == s_key_i+1).astype(int)
                        vmap_multiply = np.squeeze(vmap_multiply)[:]
                        zeroed_envelope = vmap_multiply * ae_sens
                        
                        # Find the indices that are non-zero
                        idx_nonzero = np.nonzero(zeroed_envelope)

                        # Make lists of separate beta event amplitude values
                        nonzero_envelope = envelope_cutted_list(
                            idx_nonzero, zeroed_envelope
                        )

                        # Calculate the quantiles
                        # The quantiles should be calculated to each event separately
                        cut_values = [np.quantile(non_env, cut_percentage/100) for non_env in nonzero_envelope]
                        
                        # Mean of the sections
                        mean_of_events = [
                            np.mean(seq[seq >= cut_values[i]])
                            for i, seq in enumerate(nonzero_envelope)
                        ]

                        if mean_of_events: #if the array is not empty
                            dict_append_mean[sens_id] = np.mean(mean_of_events)
                            dict_append_median[sens_id] = np.median(mean_of_events)
                            dict_append_std[sens_id] = np.std(mean_of_events)
                            dict_append_max[sens_id] = np.max(mean_of_events)

                            # Sorting for robust maximum
                            sort_index = np.argsort(mean_of_events)
                            sorted_non = np.array(mean_of_events)[sort_index]
                            rslt = sorted_non[-int(np.ceil(0.05 * len(sorted_non))) :]

                            dict_append_robust_max[sens_id] = np.mean(rslt)
                    
                    # Add then to the saving dictionary
                    save_dict.loc[len(save_dict)] = dict_append_mean
                    save_dict.loc[len(save_dict)] = dict_append_median
                    save_dict.loc[len(save_dict)] = dict_append_std
                    save_dict.loc[len(save_dict)] = dict_append_max
                    save_dict.loc[len(save_dict)] = dict_append_robust_max
                    save_dict.loc[len(save_dict)] = dict_append_whole_ts
                    save_dict.loc[len(save_dict)] = dict_append_whole_ts_low
                    save_dict.loc[len(save_dict)] = dict_append_whole_ts_high

# Save the event amplitude csv
save_dict.to_csv(fname.feature_csv(feature = 'EA',job_id=group_id, task=task))
