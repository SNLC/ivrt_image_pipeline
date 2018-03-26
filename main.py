#!/usr/bin/env python

import numpy as np
import os
import sklearn.neighbors
import pandas as pd
import sys
import time
import networkx as nx
import matplotlib.pyplot as plt


def csv_to_pandas(csv_fpath_list):
    """
    Input: List of N csv file paths
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


def test_df(df):
    valid = True
    for item in ['x', 'y']:
        if item not in df.columns:
            valid = False
            example_df = pd.DataFrame(np.array([[0, 1],
                                                [1, 2]]),
                                      columns=('x', 'y'))
            print('Warning: Column {} not in DataFrame'.format(item),
                  'Please ensure original XY csv for each',
                  'microscope channel is formatted in the',
                  'following way\n', example_df)

    for item in ['rabies', 'helper']:
        if item not in df.columns:
            valid = False
            print('Warning: Missing s and '.format(item))
    if not valid:
        sys.exit()


def multichannel_panda_to_alignment(df, radius):
    """
    Input 1: DataFrame with cell XY coordinates in first two columns
    and one hot feature vector in subsequent columnns.
    Input 2: max allowed pixel distance between two points in XY space
    for those two points to correspond to the same object
    Output: DataFrame with XY coordinates grouped if within the radius
    aligns microscope images into 'objects' possesing a one-hot vector
    which encodes if the object had the feature represented in each 
    microscope channel
    """
    print('Channel Alignment Started')
    start = time.time()
    coordinates = df.loc[:, ('x', 'y')].values
    tree = sklearn.neighbors.KDTree(coordinates, leaf_size=40)
    alignment_array = tree.query_radius(coordinates, r=radius)
    # need to fix column list to make it automatic
    column_list = ['x', 'y', 'helper', 'rabies', 'vglut']
    aligned_df = pd.DataFrame(index=df.index, columns=column_list)
    aligned_df.fillna(value=0, inplace=True)

    for idx, alignment in enumerate(alignment_array):
        xy = np.mean(coordinates[alignment],
                     axis=0,
                     dtype=np.float32)
        aligned_df.loc[idx, ('x', 'y')] = xy
        alignment_features = df.loc[alignment, 'channel_name']
        aligned_df.loc[idx, alignment_features] = 1

    pruned_aligned_df = aligned_df.drop_duplicates(keep='first')

    end = time.time()
    print("Alignment completed in {0:.2f} seconds".format(end-start))
    return pruned_aligned_df


def alignment_df_to_region_df(df, region_radius, presynaptic_radius):
    starter_cell_df = df[((df['helper'] == 1) & (df['rabies'] == 1))]
    starter_cell_xy = starter_cell_df.loc[:, ['x', 'y']]
    starter_cell_idx_array = starter_cell_df.index.values

    starter_tree = sklearn.neighbors.KDTree(starter_cell_xy,
                                            leaf_size=40)
    region_array = starter_tree.query_radius(starter_cell_xy,
                                             r=region_radius)
    region_set = {tuple(region) for region in region_array}

    region_graph = nx.Graph()
    for region in region_set:
        region_graph.add_cycle(list(region))

    joined_region_list = nx.connected_components(region_graph)

    df['starter'] = 0
    df['in region'] = 0
    df['Region ID'] = np.nan

    for region_id, region in enumerate(joined_region_list):
        indices = starter_cell_idx_array[np.array(list(region), dtype=np.int)]
        df.loc[indices, ['starter', 'in region', 'Region ID']] = [1, 1, region_id + 1]
    df.to_csv('~/Desktop/out.csv')

    unlabeled_cell_bool = df['Region ID'].isnull()
    unlabeled_cells_xy = df.loc[unlabeled_cell_bool, ['x', 'y']]
    distance_array, neighbor_array = starter_tree.query(unlabeled_cells_xy,
                                                        k=1,
                                                        dualtree=True)
    df.loc[unlabeled_cell_bool, 'distance_to_closest_starter'] = distance_array
    df.loc[unlabeled_cell_bool, 'closest_starter'] = starter_cell_idx_array[neighbor_array]
    df.loc[unlabeled_cell_bool, 'Closest Region ID'] = df.loc[df.loc[unlabeled_cell_bool, 'closest_starter'], 'Region ID'].values

    df['presynaptic'] = df['rabies'] & ~df['helper']

    proximity_bool = df['distance_to_closest_starter'] < presynaptic_radius
    proximity_mask = (unlabeled_cell_bool & proximity_bool)
    df.loc[proximity_mask, 'Region ID'] = df.loc[proximity_mask, 'Closest Region ID']
    df.to_csv('~/Desktop/region.csv')
    sys.exit()
    return df

def main():
    base_dir = ('/Volumes/My Book/rabies_tracing_images/'
                'pv_cre_starter_cells/XH_12_07_17/'
                'vglut_presynaptic_647/01_counts/')
    #base_dir = ('/Users/michael/Desktop/test/Real/')
    csv_list = ['helper.csv', 'rabies.csv', 'vglut.csv']
    csv_paths = [os.path.join(base_dir, filename)
                 for filename in csv_list]

    multichannel_df = csv_to_pandas(csv_paths)

    alignment_radius = 4.0
    alignment_df = multichannel_panda_to_alignment(multichannel_df,
                                                   alignment_radius)
    max_presynaptic_distance = 300
    region_radius = 2.0 * max_presynaptic_distance

    region_df = alignment_df_to_region_df(alignment_df,
                                          region_radius,
                                          max_presynaptic_distance)

    aligned_out_fpath = os.path.join(base_dir, 'alignment.csv')
    region_out_fpath = os.path.join(base_dir, 'region.csv')
    quant_out_fpath = os.path.join(base_dir, 'quantification.csv')


if __name__ == "__main__":
    main()
