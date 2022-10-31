import pyarrow.dataset as ds
from tqdm import tqdm


def stream(dataset_lis, file_format="parquet", batch_size=20000, columns=["canonical_ID", "enumerated_smiles", "achiral_fp"]):
    dataset = ds.dataset(dataset_lis, format=file_format)
    num_files = len(dataset.files)
    return tqdm(dataset.to_batches(columns=["canonical_ID", "enumerated_smiles", "achiral_fp"]), total=num_files)

