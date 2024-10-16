#!/usr/bin/env python
# coding: utf-8


import pandas as pd
from itertools import combinations
import sys
import argparse


def merge_bartender_output(input_barcode_address,input_cluster_address):
    temp_df1 = pd.read_csv(input_barcode_address)
    temp_df2 = pd.read_csv(input_cluster_address)
    temp_merge=pd.merge(temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
         left_index=False, right_index=False, sort=True, copy=True, indicator=False,
         validate=None)
    return(temp_merge)

def unique_read_to_cluster_dic(input_df): # the input is a merged df from two barcode.csv and cluster.csv
    temp_UniqueRead_To_ClusterSeq_dic = {}
    for x,y in zip(input_df['Unique.reads'].to_list(),input_df['Center'].to_list()):
        temp_UniqueRead_To_ClusterSeq_dic[x] = y
    return(temp_UniqueRead_To_ClusterSeq_dic)

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def All_Pairwise_Hamming_Distannce_from_df(input_df):
    # generate string pair:
    temp1 = list(combinations(input_df.gRNA, 2))
    temp_dict = {}
    for x in temp1:
        temp_dict[x] = hamming_distance(x[0],x[1])
    return (temp_dict)

def Check_sgRNA_Library_Distance(input_sgRNA_df, input_minimal_distance):
    # this module is too check if any of the two sgRNAs are within input_minimal_distance
    # for example, if we want to torelate sgRNAs with 2 PCR errors then the minimal distance should be 4
    df_groups = input_sgRNA_df.groupby("gRNA_length") # first group them based on length
    top_dict = {} # each key is the sgRNA length, value is a sub_dict
    # for each sub_dict, the key is the pair of sequence, the value is the distance
    for key, sub_df in df_groups:
        temp_dict = All_Pairwise_Hamming_Distannce_from_df(sub_df)
        temp_filtered_dict = {}
        for temp_key, value in temp_dict.items():
            if value <= input_minimal_distance:
                temp_filtered_dict[temp_key] = value
        top_dict[key] = temp_filtered_dict
    return(top_dict)

def Generate_Filtered_df(input_sgRNA_bartender_address1,input_barcode_bartender_address,
                      input_sgRNA_dic):
    # input_bartender1 and input_bartender2 should have same sequence id order. The read should be matched.
    temp_dic = {} # key is a (bacode1_cluster_seq, barcode2_cluster_seq), value is the count
    temp_sampleID = input_sgRNA_bartender_address1.split('/')[-2]
    with open(input_sgRNA_bartender_address1, 'r') as handler1, open(input_barcode_bartender_address, 'r') as handler2:
        temp1 = handler1.readline().strip()            
        temp2 = handler2.readline().strip()
        # the following lists store reads with matched sgRNA sequence
        temp_sgRNA_list, temp_sgRNA_center_list, temp_barcode_list, temp_readsID_list = [], [], [], []
        while bool(temp1)&bool(temp2):
            temp_gRNA = temp1.split(',')[0] # sgRNA 
            if temp_gRNA in input_sgRNA_dic.keys():
                temp_sgRNA_center = input_sgRNA_dic.get(temp_gRNA)
                temp_readsID = temp1.split(',')[1].split(' ')[0]
                temp_barcode = temp2.split(',')[0] # barcode
                temp_sgRNA_list.append(temp_gRNA)
                temp_sgRNA_center_list.append(temp_sgRNA_center)
                temp_barcode_list.append(temp_barcode)
                temp_readsID_list.append(temp_readsID)
            temp1 = handler1.readline().strip()            
            temp2 = handler2.readline().strip()
    temp_output_df = pd.DataFrame({'gRNA': temp_sgRNA_list, 'gRNA_center': temp_sgRNA_center_list,
                                   'Clonal_barcode':temp_barcode_list, 'Read_ID': temp_readsID_list,
                                   'Sample_ID':temp_sampleID})
    return (temp_output_df)
        
def main():
    parser = argparse.ArgumentParser(description='A function to map sgRNA clustering result to reference and seperate clonal barcode data based on sgRNA')
    parser.add_argument("--a1", required=True, help="This is the input file of bacode.csv for sgRNA")
    parser.add_argument("--a2", required=True, help="This is the input file of cluster.csv for sgRNA")
    parser.add_argument("--a3", required=True, help="This is the input file for sgRNA reference")
    parser.add_argument("--a4", required=False, help="This is the expected minimal distance of two sgRNAs")
    parser.add_argument("--a5", required=True, help="This is the input file of bartender of sgRNA")
    parser.add_argument("--a6", required=True, help="This is the input file of bartender of clonal barcode")
    parser.add_argument("--o", required=True, help="This is the prefix of output file")
    args = parser.parse_args()
    temp_output_prefix = args.o
    df1 = merge_bartender_output(args.a1,args.a2) # clustered sgRNA dataframe
    gRNA_df = pd.read_csv(args.a3)
    gRNA_df['gRNA_length'] = gRNA_df['gRNA'].apply(lambda x: len(x))
    gRNA_df = gRNA_df.drop_duplicates(keep='first') # drop duplicates
    
    if args.a4 is None: # Not checking distance
        sample_to_exclude = []
        print(f"The sgRNA reference distance is not verified")
    else:
        # check the distance of reference sgRNAs
        sgRNA_distance_dict = Check_sgRNA_Library_Distance(gRNA_df,int(args.a4)) 
        for key, value in sgRNA_distance_dict.items():
            if value != {}:
                print('sgRNA distance error\n')
                sys.exit("Error message")
        print('sgRNA are good')

    # I only look at sgRNA that matched to the reference
    df1_final = df1.copy()
    print('----------------')
    df1_final['Center'] = df1_final.apply(lambda x: x['Unique.reads']if (x['Unique.reads'] in gRNA_df.gRNA.values) else x['Center'], axis=1)
    df1_final = df1_final[df1_final.Center.isin(gRNA_df.gRNA)] # only look gRNA center that matches to the reference
    gRNA_ref_dict = unique_read_to_cluster_dic(df1_final)
    Final_df = Generate_Filtered_df(args.a5, args.a6, gRNA_ref_dict) # only saving reads with gRNA center that matches to the reference
    
    # gRNA, gRNA_center, Clonal_barcode, Read_ID, Sample_ID
    temp_sample_ID = Final_df.Sample_ID[0]
    count1 = df1.Frequency.sum()
    count2 = df1_final.Frequency.sum()
    print(f"Sample {temp_sample_ID:s} has totally {count1:d} reads for sgRNA. Among them {count2:d} match to reference sgRNA, which is about {count2/count1:.3f}")
    
    # Output intermediate files
    Final_df.to_csv(temp_output_prefix + 'Intermediate_df.csv', index=False)
    sgRNA_groups = Final_df.groupby("gRNA_center")
    # output the bartender file address
    temp_address_file = temp_output_prefix+'Bartender_input_address'
    file_a = open(temp_address_file, 'w')
    for groups in sgRNA_groups:
        g, value = groups # g=gRNA
        temp_name = temp_output_prefix + 'Clonal_barcode/' + g + '.bartender'
        file_a.write(temp_name + '\n')
        value[['Clonal_barcode', 'Read_ID']].to_csv(temp_name, sep=',', header=False, index=False)
    file_a.close()
    
if __name__ == "__main__":
    main()  


