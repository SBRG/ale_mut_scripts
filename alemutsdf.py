import os
import pandas as pd
from tqdm import tqdm


def get_mut_df(CSV_file_path,
			   include_dups,
               intragenic_muts_only):

    # Step 1: Import database
    raw_db = pd.read_csv(CSV_file_path)

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


def get_all_mut_df(dir_path,
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
    return mutat_df
