
set -e

DATA_DIR=preprocessing/data2/
TFR_DIR=preprocessing/data2/TFRecords/

printf "\n\n  CREATE TFRecord FILES..."
mkdir -p $TFR_DIR
python preprocessing/PDB2TFRecord.py \
    -annot $DATA_DIR/alphafold_annot.tsv \
    -prot_list $DATA_DIR/alphafold_train.txt \
    -npz_dir $DATA_DIR/annot_pdb_chains_npz/ \
    -num_shards 30 \
    -num_threads 30 \
    -tfr_prefix $TFR_DIR/PDB_GO_train \

python preprocessing/PDB2TFRecord.py \
    -annot $DATA_DIR/alphafold_annot.tsv \
    -prot_list $DATA_DIR/alphafold_val.txt \
    -npz_dir $DATA_DIR/annot_pdb_chains_npz/ \
    -num_shards 3 \
    -num_threads 3 \
    -tfr_prefix $TFR_DIR/PDB_GO_valid \
