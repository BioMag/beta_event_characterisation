"""
Functions required for making characteristics (files in the make_characteristics folder)
"""
import numpy as np
import mne

def beta_amplitude_envelope(data, sfreq, lower_freq, upper_freq):
    """
    Amplitude envelope obtained by using Morlet wavelet decomposition.

    Args:
        data (numpy array): The input signal data.
        sfreq (int): Sampling frequency of the data.
        lower_freq (int): Lower frequency limit.
        upper_freq (int): Upper frequency limit.

    Returns:
        (numpy array): TFR transform of the data using Morlet wavelet.
    """
    # split data into consecutive epochs
    window = 10  # length of individual time windows (in seconds)
    ws = int(window * sfreq)  # number of samples per window
    overlap = (
        1 - 0
    )  # set amount of overlap for consecutive FFT windows (second number sets amount of overlap)

    # separate data into consecutive data chunks (episode-like, because spectral_connectivity expects epochs)
    array1 = list()
    start = 0
    stop = ws
    step = int(ws * overlap)
    while stop <= data.shape[1]:
        tmp = data[:, start:stop]
        start += step
        stop += step
        array1.append(tmp)

    # define frequencies of interest
    freqs = np.arange(lower_freq, upper_freq, 1.0)
    n_cycles = freqs / 2.0

    # calculate power
    power = mne.time_frequency.tfr_array_morlet(
        array1, sfreq=sfreq, freqs=freqs, n_cycles=n_cycles, output="complex", n_jobs=16
    )

    return power, freqs


def amplitude_envelope(lower_freq, higher_freq, freqs, power):
    """This function averages over powers of multiple frequencies 

    Args:
        lower_freq (int): Lower frequency limit.
        higher_freq (int): Upper frequency limit.
        freqs (numpy array): all possible frequencies
        power (numpy array): Power from the Morlet wavelet decomposition

    Returns:
        (numpy array): Amplitude envelope for each channels
    """

    b_idx = np.where(np.logical_and(freqs > lower_freq, freqs < higher_freq))
    amplitude = np.mean(np.abs(power[0][:, b_idx[0]]), axis=1)

    for k in range(1, len(power)):
        tmp = power[k][:, b_idx[0]]
        tmptmp = np.mean(np.abs(tmp), axis=1)
        amplitude = np.concatenate((amplitude, tmptmp), axis=1)

    return amplitude


def envelope_cutted_list(idx_nonzero, zeroed_envelope):
    """Segments an amplitude envelope into contiguous parts based on non-zero indices.

    Args:
        idx_nonzero (numpy array): list of indices that have nonzero values
                                    in the zeroed envelope list
        zeroed_envelope (numpy array):  Array representing the amplitude envelope
                                    values from which non-zero segments are extracted.

    Returns:
        list: List of numpy arrays representing segments of defined by non-zero indices.
    """
    fi = 0
    empty_list = []
    for ind in range(1, len(idx_nonzero[0])):
        if (idx_nonzero[0][ind] - idx_nonzero[0][ind - 1] != 1) and (
            ind != len(idx_nonzero[0])
        ):
            take_indices = idx_nonzero[0][fi:ind]
            take_envelope = zeroed_envelope[take_indices]
            empty_list.append(take_envelope)
            fi = ind

        elif ind == len(idx_nonzero[0]):
            take_indices = idx_nonzero[0][fi : ind + 1]
            take_envelope = zeroed_envelope[take_indices]
            empty_list.append(take_envelope)
    return empty_list


def cut_ones_into_list(vmap_multiply):
    """
    Takes vmap_multiply in and separates that into separate
    lists having zeros and ones.

    Args:
        vmap_multiply (numpy array): Array in which 1=state is on and 0=state is off

    Returns:
        list: lists having the sequencies of zeros and ones separately
    """

    fi = 0
    list_ones = []
    list_zeros = []
    for ind in range(1, len(vmap_multiply)):
        if vmap_multiply[ind] + vmap_multiply[ind - 1] == 1:
            take_list = vmap_multiply[fi:ind]
            if take_list[0] == 0:
                list_zeros.append(take_list)
            else:
                list_ones.append(take_list)
            fi = ind
    return list_ones, list_zeros


def AE_per_event(viterbi_path, amplitude_envelope):
    """separate viterbi path and amplitude envelope into separate events

    Args:
        viterbi_path (numpy array): viterbi path, the state sequence
        amplitude_envelope (numpy array): amplitude envelope of the signal

    Returns:
        numpy array: (sequencies of the states, signal belonging to that sequence, indexes of the sequence)
    """
    
    
    separated_signals = []
    state_seqs = []
    idx_seqs = []

    current_ids = [0]
    current_group = [amplitude_envelope[0]]
    current_value = viterbi_path[0]

    idxs = list(range(len(viterbi_path))) 

    for num, sig, id in zip(viterbi_path[1:], amplitude_envelope[1:],idxs[1:]):
        if num == current_value:
            current_group.append(sig)
            current_ids.append(id)
        else:
            separated_signals.append(np.array(current_group))
            idx_seqs.append(np.array(current_ids))
            state_seqs.append(int(current_value))
            
            current_ids = [id]
            current_group = [sig]
            current_value = num
    
    separated_signals.append(np.array(current_group))
    idx_seqs.append(np.array(current_ids))
    state_seqs.append(int(current_value))
    
    return state_seqs,separated_signals, idx_seqs



def SLT_per_event(viterbi_path):
    """separate viterbi path into events, so that it returns the state sequencies and their respective indices

    Args:
        viterbi_path (numpy array): viterbi path, the state sequence

    Returns:
        numpy array: (sequencies of the states, indexes of the sequence)
    """

    state_seqs = []
    idx_seqs = []

    current_ids = [0]
    current_value = viterbi_path[0]

    idxs = list(range(len(viterbi_path))) 

    for num, id in zip(viterbi_path[1:],idxs[1:]):
        if num == current_value:
            current_ids.append(id)
        else:
            idx_seqs.append(np.array(current_ids))
            state_seqs.append(int(current_value))
            
            current_ids = [id]
            current_value = num
    
    idx_seqs.append(np.array(current_ids))
    state_seqs.append(int(current_value))
    
    return state_seqs, idx_seqs