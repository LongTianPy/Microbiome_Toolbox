#!/usr/bin/python
"""
"""

# IMPORT
import pandas as pd
import numpy as np
import sys

# VARIABLES


# FUNCTIONS
def convert(observed_otu_file,mapping_file,by):
    """
    This function convert the observed otus file from QIIME to a format that R can take and draw rarefaction curve
    :param observed_otu_file:
    :return:
    """
    df = pd.read_table(observed_otu_file,sep="\t",header=0,index_col=0)
    sample_id = df.columns[2:]
    depths = list((set(list(df['sequences per sample']))))
    depths.sort()
    Sample = []
    Depth = []
    Mean = []
    Std = []
    mapping = pd.read_table(mapping_file,sep="\t",header=0,index_col=0)
    mapping["SampleID"] = list(mapping.index)
    group_by = {}
    by_group = {}
    for i in mapping.index:
        group_by[i] = mapping.loc[i,by]
        if mapping.loc[i,by] not in by_group:
            by_group[mapping.loc[i,by]]=[i]
        else:
            by_group[mapping.loc[i,by]].append(i)
    for depth in depths:
        sub_df_by_depth =  df[df["sequences per sample"]==depth].iloc[:,2:]
        for each_group in by_group:
            sub_sub = sub_df_by_depth.loc[:,by_group[each_group]]
            Sample.append(each_group)
            Depth.append(depth)
            sample_mean = np.mean(sub_sub.values)
            sample_std = np.std(sub_sub.values)
            Mean.append(sample_mean)
            Std.append(sample_std)
    ## Simple version
    # for depth in depths:
    #     sub_df = df[df["sequences per sample"]==depth]
    #     for each_sample in sample_id:
    #         Sample.append(each_sample)
    #         Depth.append(depth)
    #         sample_mean = np.mean(sub_df[each_sample])
    #         sample_std = np.std(sub_df[each_sample])
    #         Mean.append(sample_mean)
    #         Std.append(sample_std)
    new_df = pd.DataFrame()
    new_df["Sample"] = Sample
    new_df["Depth"] = Depth
    new_df["Mean"] = Mean
    new_df["Std"] = Std
    new_df.to_csv("transformed_observed_otus.txt",sep="\t")

# MAIN
if __name__ == '__main__':
    inputfile = sys.argv[1]
    mappingfile = sys.argv[2]
    by = sys.argv[3]
    convert(inputfile,mappingfile,by)