import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
import glob
import numpy as np
import pyarrow.feather as feather
# import dataframe_image as dfi
# from rdkit.Chem import PandasTools
import sys
import molbeam
import similarity

dataset_dir = sys.argv[1]
output_dir = sys.argv[2]
# out_prefix = sys.argv[3]

def process_batch_fp(query_names, query_matrix, mol_batch):
    fingerprints = mol_batch.column('achiral_fp')
    fp_matrix_in = similarity.build_fingerprint_matrix(fingerprints)
    fp_distance = similarity.fast_jaccard(fp_matrix_in, query_matrix)
    # fp_distance_matrix = np.asmatrix(fp_distance)
    # print(type(fp_distance))

    # # if elements are not below threshold don't keep them
    # matches = fp_distance_matrix <= threshold
    # result = None
    # if matches.any():
    #     result = pd.DataFrame(fp_distance, columns=['PBB3', 'PM_PBB3', 'C5_05', 'chet_hit1', 'chet_hit2', 'chet_hit3', 'chet_hit4'])
    #     result.insert(0, 'canonical_id', mol_batch.column('canonical_ID'))
    #     result.insert(1, 'std_smiles', mol_batch.column('enumerated_smiles'))
    #     print('found matches')

    return fp_distance,  mol_batch.column('enumerated_smiles'), mol_batch.column('canonical_ID')

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
    result_df = result_df.sort_values(by=['name', 'name', 'name', 'name', 'name','name', 'name'], ascending=True)
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
    # minimum quatile to be considered a match
    quantile = 0.00005
    results = []
    for mol_batch in molbeam.stream(lis, file_format="parquet", batch_size=20000, columns=columns):
        result = process_batch_fp(query_names, query_matrix, mol_batch)
        results.append(result)


    scores_lis = [results[_][0] for _ in range(len(results))]
    scores_arr = np.vstack(scores_lis)
    smiles_lis_pa = [results[_][1] for _ in range(len(results))]
    canonical_id_lis_pa = [results[_][2] for _ in range(len(results))]
    smiles_arr = pa.concat_arrays(smiles_lis_pa)
    canonical_id_arr = pa.concat_arrays(canonical_id_lis_pa)
    d = {}
    for i,query in enumerate(query_names):
        query_col = f'{query}_score'
        scores_col = scores_arr[:, i]
        query_quantile = np.quantile(scores_col, quantile)
        d[query_col] = {
            'scores':scores_col,
            'quantile':query_quantile 
        }
    # pbb3_scores = pa.array(scores_arr[:, 0])
    # pm_pbb3_scores = pa.array(scores_arr[:, 1])
    # c505_scores = pa.array(scores_arr[:, 2])
    data = [
        canonical_id_arr,
        smiles_arr,
    ]
    col_names = [
        'canonical_id',
        'enumerated_smiles',
    ]
    for i,col_name in enumerate(d):
        arr = d[col_name]['scores']
        data.append(arr)
        col_names.append(col_name)
    table = pa.Table.from_arrays(data, names=col_names)
    # feather.write_feather(
    #     table, 
    #     "/scratch/grahamth/test.feather", 
    #     compression="uncompressed", 
    #     chunksize=1000_000
    # )

    # source = pa.memory_map("/scratch/grahamth/test.feather", 'r')
    # reader = pa.RecordBatchFileReader(source)
    # table = reader.read_all()
    # df = 
    # table.to_pandas()
    # export_results(results, query, threshold=threshold)
    rel = duckdb.from_arrow_table(table)
    sql_lis = []
    for i,col_name in enumerate(d):
        quantile = d[col_name]['quantile']
        sql = '''
            SELECT *
            FROM arrow
            WHERE {col_name}
            < {quantile}
        '''.format(col_name=col_name, quantile=quantile)
        sql_lis.append(sql)


    df_lis = []
    for sql in sql_lis:
        df = rel.query("arrow", sql).df()
        # PandasTools.AddMoleculeColumnToFrame(df, smilesCol='enumerated_smiles', molCol='hit')
        # PandasTools.ChangeMoleculeRendering(df, renderer="SVG")
        df_lis.append(df)


    bigdata = pd.concat(df_lis, ignore_index=True, sort=False)
    bigdata.reset_index()
    print(len(bigdata))
    table = pa.Table.from_pandas(bigdata)
    pq.write_table(table, f'{output_dir}/query_out.parquet')
    print(f'Wrote results to: {output_dir}/query_out.parquet')

if __name__ == '__main__':
    main()
