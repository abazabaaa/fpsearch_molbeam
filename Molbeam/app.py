import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
import glob
import numpy as np
# import dataframe_image as dfi
# from rdkit.Chem import PandasTools
import sys
import molbeam
import similarity

dataset_dir = sys.argv[1]
output_dir = sys.argv[2]
# out_prefix = sys.argv[3]

def process_batch_fp(query_names, query_matrix, mol_batch, threshold=.99):
    fingerprints = mol_batch.column('achiral_fp')
    fp_matrix_in = similarity.build_fingerprint_matrix(fingerprints)
    fp_distance = similarity.fast_jaccard(fp_matrix_in, query_matrix)
    fp_distance_matrix = np.asmatrix(fp_distance)
    # print(type(fp_distance))

    # if elements are not below threshold don't keep them
    matches = fp_distance_matrix <= threshold
    result = None
    if matches.any():
        result = pd.DataFrame(fp_distance, columns=['xxx', 'xxx', 'xxx', 'xxx', 'xxxx', 'xxx', 'xxx'])
        result.insert(0, 'canonical_id', mol_batch.column('canonical_ID'))
        result.insert(1, 'std_smiles', mol_batch.column('enumerated_smiles'))
        print('found matches')

    return result

# def process_batch_fastrocs(query_names, query_matrix, mol_batch, threshold=.99):
    # fingerprints = mol_batch.column('achiral_fp')
    # fp_matrix_in = similarity.build_fingerprint_matrix(fingerprints)
    # fp_distance = similarity.fast_jaccard(fp_matrix_in, query_matrix)

    # # if elements are not below threshold don't keep them
    # matches = fp_distance <= threshold
    # result = None
    # if matches.any():
    #     result = pd.DataFrame(fp_distance, columns=['PBB3', 'PM_PBB3', 'C5_05'])
    #     result.insert(0, 'canonical_id', mol_batch.column('canonical_ID'))
    #     result.insert(1, 'std_smiles', mol_batch.column('enumerated_smiles'))
    #     print('found matches')

    # return result


def export_results(result_list, query, threshold):
    result_df = pd.concat(result_list)
    result_df = result_df.sort_values(by=['xxx', 'xxx', 'xxx', 'xx', 'xx','xx', 'xxx'], ascending=True)
    print(result_df.head(10))
    table = pa.Table.from_pandas(result_df)
    pq.write_table(table, f'{output_dir}/query_out.parquet')
    print(f'Wrote results to: {output_dir}/query_out.parquet')

    # con = duckdb.connect(database=':memory:', read_only=False)
    # sql = '''
    #     SELECT
    #         *
    #     FROM parquet_scan('{output_dir}/query_out.parquet')
    #     WHERE PBB3 < ?
    # '''.format(output_dir=output_dir)
    # top_results = con.execute(sql, [threshold]).fetchdf()

    # # PandasTools.AddMoleculeColumnToFrame(top_results, smilesCol='std_smiles', molCol='self')
    # final_df = top_results.sort_values(by=['PBB3'], ascending=True)
    # final_df = final_df.reset_index(drop=True)
    # print(final_df.head(10))
    # table_hits = pa.Table.from_pandas(final_df)
    # pq.write_table(table_hits, f'{output_dir}/query_out_hits.parquet')
    # print('Wrote results to: query_out_hits.parquet')
    # file_name = 'search_results.png'
    # dfi.export(final_df, file_name)


def main():
    query = [
        ('name', 'smiles'),
        ('name', 'smiles'),
        ('name', 'smiles'),
        ('name', 'smiles'),
        ('name', 'smiles'),
        ('name', 'smiles'),
        ('name', 'smiles')
    ]
    query_names = [q[0] for q in query]
    query_matrix = similarity.format_query(query)
    columns = ["canonical_ID", "enumerated_smiles", "achiral_fp"]
    results = []
    # print('Searching enamine database of 3.8M molecules...')

    dataset = dataset_dir
    lis = glob.glob(f'{dataset}/*.parquet')
    # minimum jaccard distance to be considered a match
    threshold = 0.75
    for mol_batch in molbeam.stream(lis, file_format="parquet", batch_size=20000, columns=columns):

        result = process_batch_fp(query_names, query_matrix, mol_batch, threshold)
        if result is not None:
            results.append(result)



    export_results(results, query, threshold=threshold)


if __name__ == '__main__':
    main()