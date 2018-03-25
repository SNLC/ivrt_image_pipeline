#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import sklearn.neighbors
import pandas as pd

def csv_to_pandas(csv_fpath_list):
    """Input: List of N csv file paths and max allowed pixel distance 
    between two points in XY space for those two points to correspond
    to the same object
    Output: Pandas DataFrame where each row corresponds to a single cell,
    the first two columns correspond to mean X and Y coordinates between
    aligned points. The remaining columns are named after the CSV files
    and take binary truth values if that point has the 'feature' corresponding
    to that microscope channel"""
    
    frames = []
    for csv_fpath in csv_fpath_list:
        df = pd.read_csv(csv_fpath) 
        channel_name = os.path.splitext(os.path.basename(csv_fpath))[0]
        df['channel_name'] = channel_name 
        frames.append(df)

    df = pd.concat(frames)
    df.reset_index(drop=True, inplace=True)
    return df 

def multichannel_panda_to_alignment(df, radius):
    coordinates = df.loc[:,('x','y')].values
    tree = sklearn.neighbors.KDTree(coordinates, leaf_size=40)
    alignment_array = tree.query_radius(coordinates, r=radius)
    #need to fix column list to make it automatic
    column_list = ['x','y','helper','rabies', 'vip']
    aligned_df = pd.DataFrame(index=df.index, columns=column_list)
    
    for idx, alignment in enumerate(alignment_array):
        xy = np.mean(coordinates[alignment],
                axis=0,
                dtype=np.float32)
        aligned_df.loc[idx, ('x', 'y')] = xy 
        alignment_features = df.loc[alignment, 'channel_name']
        aligned_df.loc[idx, alignment_features] = 1

    pruned_df = aligned_df.drop_duplicates(keep='first')
    print(pruned_df)
    return pruned_df

def main():
    base_dir = ('/Volumes/My Book/rabies_tracing_images/'
               'pv_cre_starter_cells/XH_12_07_17/'
               'vip_presynaptic_647/01_counts/')
    csv_list = ['helper.csv', 'rabies.csv', 'vip.csv']
    csv_paths = [base_dir + filename for filename in csv_list]

    multichannel_df = csv_to_pandas(csv_paths) 

    alignment_radius = 4.0
    multichannel_panda_to_alignment(multichannel_df, alignment_radius)
        
    max_presynaptic_distance = 300
    region_radius = 2.0 * max_presynaptic_distance
    aligned_out_fpath = os.path.join(base_dir, 'alignment.csv')
    quant_out_fpath = os.path.join(base_dir, 'alignment.csv')
main()
