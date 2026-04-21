"""
Functions required for preprocessing (files in the process_raw_data folder)
"""
import mne


# Make run ICA function
def run_ica(method, raw_or_epochs, n_components = 15):
    """_summary_

    Args:
        method (string): method used to calculate the ICA decomposition.
        raw_or_epochs (Epochs or Raw): The data to be .
        n_components (int, optional): Number of components who wants to make. Defaults to 15.

    Returns:
        (ICA): ICA solution
    """
    ica = mne.preprocessing.ICA(n_components=n_components, method=method, 
              random_state=97, max_iter=400)
    ica.fit(raw_or_epochs, picks = ['meg'])
    
    return ica


def reject_bad_segments(raw):
    """
    Function to reject bad data segments from raw data

    Args:
        raw (Raw): raw data

    Returns:
        raw_out (Raw): data without bad segments
    """

    if len(raw.annotations) >= 1:
        fs = raw.first_samp
        sf = raw.info['sfreq']

        raw_mins = [0]
        raw_maxs = [(raw.n_times/sf)-(1/sf)]

        for j in range(0, len(raw.annotations)):
            
            tmax = raw.annotations.onset[j] # Onset to the max list
            tmin = raw.annotations.onset[j] + raw.annotations.duration[j] # End to the min list
            raw_mins.append(tmin - (fs/sf))
            raw_maxs.append(tmax - (fs/sf))

        # Sort the start and end lists
        raw_mins.sort()
        raw_maxs.sort()

        # Zip the lists to get tmins and maxs to crop function
        crop_sections = list(zip(raw_mins,raw_maxs))
        
        # Empty list for raw segments
        raw_segs = []

        for crop in crop_sections:
            raw_segs.append( raw.copy().crop(tmin=crop[0],tmax=crop[1]))
        
        raw_out = mne.concatenate_raws(raw_segs)
        return raw_out

    else:
        raw_out = raw
        return raw_out
