import os
import pandas as pd
from tqdm import tqdm


def get_mut_dataframe(CSV_file_path,
                      include_dups,
                      intragenic_muts_only):
    
    # Step 1: Import database
    raw_db = pd.read_csv(CSV_file_path)
    if 'Function' in raw_db.columns:
        raw_db = raw_db.drop('Function', axis=1)
    if 'Product' in raw_db.columns:
        raw_db = raw_db.drop('Product', axis=1)
    if 'GO Process' in raw_db.columns:
        raw_db = raw_db.drop('GO Process', axis=1)
    if 'GO Component' in raw_db.columns:
        raw_db = raw_db.drop('GO Component', axis=1)
    
    # Step 2: Separate columns based on usage       
    keep_cols = ['Position','Mutation Type','Sequence Change','Details','Gene']
    mut_cols = sorted(list(set(raw_db.columns) - set(keep_cols)))

    # Step 3: Shift mutation column names into row identifiers
    csv_file_mutat_df = pd.DataFrame()
    for col in tqdm(mut_cols):
        df = raw_db[raw_db[col].notnull()][keep_cols]
        df['exp'] = '_'.join(col.split(' ')[:-4])
        df['ale'] = int(col.split(' ')[-4][1:])
        df['flask'] = int(col.split(' ')[-3][1:])
        df['isolate'] = int(col.split(' ')[-2][1:])
        df['tech_rep'] = int(col.split(' ')[-1][1:])
        df['presence'] = raw_db[raw_db[col].notnull()][col]
        csv_file_mutat_df = pd.concat([csv_file_mutat_df,df])

    csv_file_mutat_df = csv_file_mutat_df[['exp','ale','flask','isolate','tech_rep','presence'] + keep_cols]
    csv_file_mutat_df = csv_file_mutat_df.fillna('')
    
    # Remove mutation entries with empty gene since they will screw up mutat_df.groupby(['Gene', ...])
    csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Gene'] != '']
    
    # Remove weird characters between gene names in multiple gene annotation.
    csv_file_mutat_df['Gene'] = csv_file_mutat_df['Gene'].str.replace("  ", " ")
    
    if not include_dups:
        csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Details'] != 'Duplication']
    
    if intragenic_muts_only:
        csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Gene'].str.contains(',') == False]
    
    return csv_file_mutat_df


def get_all_sample_mut_df(dir_path,
                          include_dups,
                          intragenic_muts_only):
    mutat_df = pd.DataFrame()
    mutat_df_list = []
    file_path_list = []
    for file_name in os.listdir(dir_path):
        print(file_name)
        file_path = dir_path+'/'+file_name
        file_path_list.append(file_path)        
        mutat_df_list.append(get_mut_dataframe(file_path, include_dups, intragenic_muts_only))
        
    mutat_df = pd.concat(mutat_df_list)
#    for file_path in file_path_list: print(file_path)
    return mutat_df


def _get_exp_ale_set(mut_df):
    exp_ale_df = mut_df.copy()
    exp_ale_df["exp ale"] = exp_ale_df["exp"] + ' ' + exp_ale_df["ale"].map(str)
    exp_ale_set = set(exp_ale_df['exp ale'].tolist())
    return exp_ale_set


def get_gene_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Sequence Change",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] = 1
    
    return mut_mat_df


def get_gene_mut_count_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] += 1
    
    return mut_mat_df


def get_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"] + ' ' + mut_df["ale"].map(str)
    mut_df["gene seq change"] = mut_df["Gene"] + ' ' + mut_df["Sequence Change"]
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)  # Get unique set of exp+ALE#

    mut_df_index_set = set(mut_df["gene seq change"].tolist())  # Get unique set of gene+(seq change) names.
    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for mut, all_mut_df in mut_df.groupby("gene seq change"):
        for index, mut_df_row in all_mut_df.iterrows():
            mut_mat_df.loc[mut, mut_df_row["exp ale"]] = 1
    
    return mut_mat_df


def get_enrichment_muts(mut_df):
    trunc_mut_df = mut_df.copy()

    # If we are going to keep Duplications, though we want to remove the '[' and ']' from their gene annotations.
#     trunc_mut_df["Gene"] = trunc_mut_df["Gene"].map(lambda x: x.lstrip('[').rstrip(']'))

    # Removing duplications
    trunc_mut_df = trunc_mut_df[trunc_mut_df["Details"] != "Duplication"]

    # Removing unused columns
    del trunc_mut_df["tech_rep"]
    del trunc_mut_df["isolate"]
    # Could have the same mutation, but with a different presence due to differences between clonal and population reseq'ing
    del trunc_mut_df["presence"]
    trunc_mut_df = trunc_mut_df.drop_duplicates()

    enrichment_mut_df = pd.DataFrame()
    for gene_mut_groupby in trunc_mut_df.groupby(["exp", "Gene"]):
        gene_mut_df = gene_mut_groupby[1]
        mutation_count = gene_mut_df.shape[0]
        if mutation_count > 1:
            enrichment_mut_df = enrichment_mut_df.append(gene_mut_df)

    return enrichment_mut_df