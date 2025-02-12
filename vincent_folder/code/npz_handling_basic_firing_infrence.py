#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:15:55 2025

@author: vincentcalia-bogan
"""

# Guidance to handling Vincent's spike-train .npz data

# function that extracts data from my npz file-- this is constructed as a generator due to RAM constraints
## MAKING A GENERATOR FOR THE NPZ FILES FOR FURTHER OPERATIONS ## 
# TODO: see if smarter to revise from generator to a polars dataframe 
# note: this is very old code...
import os
import numpy as np  
def extract_from_npz(save_path=None):
    if save_path is None:
        save_path = input("Enter the save path: ")

    # Check if save_path exists
    if not os.path.exists(save_path):
        print(f"Error: Save path '{save_path}' does not exist.")
        return

    # Check if save_path is a directory
    if not os.path.isdir(save_path):
        print(f"Error: '{save_path}' is not a directory.")
        return

    # Iterate over .npz files in save_path
    npz_files = [f for f in os.listdir(save_path) if f.endswith('.npz')]
    if not npz_files:
        print(f"No .npz files found in '{save_path}'.")
        return
    npz_files = os.listdir(save_path)
    for npz_file in npz_files:
        if npz_file.endswith('.npz'): # filter out other potential files 
            file_path = os.path.join(save_path, npz_file)
            npz_data = np.load(file_path) # load data, create a dictionary 
            for index, key in enumerate(npz_data.keys()): 
                if key.startswith('spike_array '): # indexing data correctly using a dictionary; the key command
                    spike_array = npz_data[key] # create a generator for a given spike array
                    yield spike_array, npz_file, index, key # decide if index and key are needed
                    
# Example: extracting data-- this needs to be a directory and not the file itself:  
save_path = '/Users/vincentcalia-bogan/Desktop/1BRANDEIS MAJOR STUFF/Senior Year/NBIO 207A/vincent_folder/data'
for data in extract_from_npz(save_path): # save path is the dir that you've got the .npz file in
    if isinstance(data, tuple):
        spike_array, dataset_num, index, key = data
        print(f'Dataset number: {dataset_num}')
        num_tastes, num_trials, num_neurons, num_time = spike_array.shape # construction of npz file, which is tastes, trials, neurons, time 

    for taste_idx in range(num_tastes): 
        for trial_idx in range(num_trials): 
            trial_spike_array = spike_array[taste_idx, trial_idx, :, :]
        
# example usage -- inferring firing rate via rolling window
# note this is NOT the most efficient way to do this; it is however a bit more intuitive 
# I have found that a window length of 250 ms and a step size of 25 ms works well 
# this particular function wants data constructed according to the type of trial_spike_array
def calc_fr_rr(trial_spike_array, window_length, step_size):
    num_bins = max((trial_spike_array.shape[1] - window_length) // step_size + 1, 1)
    nrn_num, time_num = trial_spike_array.shape
    firing_rate = np.zeros((nrn_num, num_bins)) # 2x2 array 
    for bin_ini in range(num_bins):
        start_bin = bin_ini * step_size
        end_bin = min(start_bin + window_length, time_num)
        num_spikes = np.sum(trial_spike_array[:, start_bin:end_bin], axis=-1)
        firing_rate[:, bin_ini] = (num_spikes / (window_length / 1000))  # Convert to Hz
    return firing_rate

# calling the function: 
fr = calc_fr_rr(trial_spike_array, window_length = 250, step_size=25)
# understanding the output of the inferred firing rate: 
num_neurons = fr.shape[0] # the number of neurons in the dataset we inferred the firing rate for 
num_bins = fr.shape[1] # the number of time-bins (250 ms each) we inferred the firing rate for 
# usage varies drastically from here...


