#!/usr/bin/env python

import numpy as np
import os
import sklearn.neighbors
import pandas as pd
import sys
import time
import networkx as nx


def test_csv(csv_list):
    """
    Input: List of N csv file paths
    Output: Exits script if file paths do not contain:
    1. CSV files with 'x' and 'y' columns
    2. A 'rabies.csv' and 'helper.csv' file
    3. At least one csv file not in 2 (target_cell_type.csv
    (ex: 'vglut.csv', 'pv.csv', 'vip.csv'...)
    """

    for csv in csv_list:
        try:
            df = pd.read_csv(csv)
        except OSError:
            print("File: {} does not exist".format(csv))
            sys.exit()

        valid = all(x in df.columns for x in ['x', 'y'])
        if not valid:
            example_df = pd.DataFrame(np.array([[0, 1],
                                                [1, 2]]),
                                      columns=('x', 'y'))
            print("Column Names 'x' and 'y' not in {}".format(csv),
                  "Please ensure CSV for each channel",
                  "is formatted in the following way",
                  "(case matters)\n",
                  example_df)
            sys.exit()

    csv_file_names = [os.path.splitext(os.path.basename(x))[0]
                      for x in csv_list]

    if not all(x in csv_file_names for x in ['rabies', 'helper']):
        print("Missing 'rabies.csv' or 'helper.csv'")
        sys.exit()

    if len(csv_file_names) <= 2:
        print("Missing target cell type file",
              "(ex: vglut.csv, 'vip.csv', 'pv.csv')")
        sys.exit()


def csv_to_pandas(csv_fpath_list):
    """
    Input: List of N csv file paths
    Output: Pandas DataFrame with rows equal to sum of rows across
    CSV files. Column 1 is the X coordinate, column 2 is the
    Y coordinate, and column 3 is the CSV basename which contained
    the coordinate
    """
    frames = []
    for csv_fpath in csv_fpath_list:
        df = pd.read_csv(csv_fpath)
        channel_name = os.path.splitext(os.path.basename(csv_fpath))[0]
        df['channel_name'] = channel_name
        frames.append(df)

    df = pd.concat(frames)
    df.reset_index(drop=True, inplace=True)
    return df


def multichannel_panda_to_alignment(df, csv_list, radius):
    """
    Input 1: DataFrame with cell XY coordinates in first two columns
    and CSV name from which the row originated.
    Input 2: max allowed pixel distance between two points in XY space
    for those two points to correspond to the same object
    Output: DataFrame with XY coordinates grouped if within a given radius.
    Aligns microscope images into 'objects' possesing a one-hot vector
    which encodes if the object had the feature represented in each
    microscope channel
    """
    print('Channel Alignment Started')
    start = time.time()
    coordinates = df.loc[:, ('x', 'y')].values
    tree = sklearn.neighbors.KDTree(coordinates, leaf_size=40)
    alignment_array = tree.query_radius(coordinates, r=radius)
    column_list = ['x', 'y'] + [os.path.splitext(csv)[0] for csv in csv_list]
    aligned_df = pd.DataFrame(index=df.index, columns=column_list)
    aligned_df.fillna(value=0, inplace=True)

    for idx, alignment in enumerate(alignment_array):
        xy = np.mean(coordinates[alignment],
                     axis=0,
                     dtype=np.float32)
        aligned_df.loc[idx, ('x', 'y')] = xy
        alignment_features = df.loc[alignment, 'channel_name']
        aligned_df.loc[idx, alignment_features] = 1

    aligned_df.drop_duplicates(keep='first', inplace=True)

    end = time.time()
    print("Alignment completed in {0:.2f} seconds".format(end-start))
    return aligned_df, column_list


def alignment_df_to_region_df(df, region_radius, presynaptic_radius):
    """
    Input 1: Aligned DataFrame with feature vectors
    Input 2: Radius within which starter cells are grouped into regions
    Input 3: Radius within which presynaptic cells are assigned to a
    given starter cell region
    Output: DataFrame containing each starter cell, presynaptic cell,
    and non-presynaptic cell grouped according to closest region.
    """
    df.loc[((df['rabies'] == 1) & (df['helper'] == 0)), 'presynaptic'] = 1
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
        df.loc[indices, ['starter', 'in region']] = [1, 1]
        df.loc[indices, 'Region ID'] = region_id + 1
    # Better name for unlabeled_cell_bool = not_starter_cell
    unlabeled_cell_bool = (df['in region'] == 0)
    unlabeled_cells_xy = df.loc[unlabeled_cell_bool, ['x', 'y']]
    distance_array, neighbor_array = starter_tree.query(unlabeled_cells_xy,
                                                        k=1,
                                                        dualtree=True)

    df.loc[unlabeled_cell_bool, 'distance_to_closest_starter'] = distance_array
    # Comment this line of code could use iloc,
    # probably easier to read (iloc(df, :) index.values)
    df.loc[unlabeled_cell_bool, 'closest_starter'] = starter_cell_idx_array[neighbor_array]

    closest_starter_cell = df.loc[unlabeled_cell_bool, 'closest_starter']
    closest_region_id = df.loc[closest_starter_cell, 'Region ID'].values
    df.loc[unlabeled_cell_bool, 'Closest Region ID'] = closest_region_id

    proximity_bool = df['distance_to_closest_starter'] < presynaptic_radius
    prox_mask = (unlabeled_cell_bool & proximity_bool)
    df.loc[prox_mask, 'Region ID'] = df.loc[prox_mask, 'Closest Region ID']
    df.loc[prox_mask, 'in region'] = 1

    return df


def regional_statistics(df, columns, target_cell_type):
    presynaptic_target_label = 'presynaptic ' + target_cell_type
    target_mask = ((df[target_cell_type] == 1) & (df['rabies'] == 1))
    df.loc[target_mask, presynaptic_target_label] = 1
    # This could be more robust... should propogate column labels in a more
    # pythonic way
    columns_to_sum = columns[2:] + ['starter', 'presynaptic', presynaptic_target_label]
    summary_df = df.set_index('Region ID')
    summary_df.sort_index(inplace=True)
    summary_df = summary_df.loc[:, columns_to_sum].sum(level=0)

    presynaptic_target_presynaptic = presynaptic_target_label + '/presynaptic'
    summary_df[presynaptic_target_presynaptic] = summary_df[presynaptic_target_label]/summary_df['presynaptic']

    presynaptic_target_target = presynaptic_target_label + '/' + target_cell_type
    summary_df[presynaptic_target_target] = summary_df[presynaptic_target_label]/summary_df[target_cell_type]

    summary_df['presynaptic per starter'] = summary_df['presynaptic']/summary_df['starter']
    return summary_df


def main():
    # base_dir = ('/Volumes/My Book/rabies_tracing_images/'
    #             'pv_cre_starter_cells/XH_12_07_17/'
    #             'vglut_presynaptic_647/01_counts/')
    base_dir = '~/Desktop/test/Synthetic/'
    csv_list = ['helper.csv', 'rabies.csv', 'vglut.csv']
    csv_paths = [os.path.join(base_dir, filename)
                 for filename in csv_list]
    test_csv(csv_paths)

    multichannel_df = csv_to_pandas(csv_paths)

    alignment_radius = .1
    alignment_df, column_list = multichannel_panda_to_alignment(multichannel_df,
                                                                csv_list,
                                                                alignment_radius)
    max_presynaptic_distance = 1 
    region_radius = 2.0 * max_presynaptic_distance

    region_df = alignment_df_to_region_df(alignment_df,
                                          region_radius,
                                          max_presynaptic_distance)

    summary_df = regional_statistics(region_df, column_list, 'vglut')

    region_df.to_csv(os.path.join(base_dir, 'region.csv'))
    summary_df.to_csv(os.path.join(base_dir, 'summary.csv'))


if __name__ == "__main__":
    main()
