# !/bin/bash
# Script for training DeepFRI model

# cuda libraries to run tf2 on gpu
# module load slurm gcc cuda/10.1.105_418.39 cudnn/v7.6.2-cuda-10.1

graph_conv_dims="512 512 512"
fully_connected_dims="1024"
graph_conv_layer=GraphConv
ontology=mf
cmap_thresh=10.0
data_dir=preprocessing/data2/
cmap_data=alphafold # possible: PDB, SWISS-MODEL or MERGED
model_name=./results/DeepFRI-${cmap_data}_${graph_conv_layer}_gcd_$(echo $graph_conv_dims | tr ' ' '-')_fcd_${fully_connected_dims}_ca_${cmap_thresh}_${ontology}

if [ "$ontology" == ec ]; then
    annot=EC
    annot_fn=./preprocessing/data2/alphafold_annot.tsv
    test_list=./preprocessing/data2/alphafold_test.csv
else
    annot=GO
    annot_fn=./preprocessing/data2/alphafold_annot.tsv
    test_list=./preprocessing/data2/alphafold_test.csv
fi

if [ "$cmap_data" == MERGED ]; then
    cmap_data=*
fi

echo $annot
python train_DeepFRI.py \
    -gcd ${graph_conv_dims} \
    -fcd ${fully_connected_dims} \
    -l2 2e-5 \
    -lr 0.0002 \
    -gc ${graph_conv_layer} \
    -e 50 \
    --batch_size 32 \
    -ont ${ontology} \
    -lm trained_models/lstm_lm_tf.hdf5 \
    --cmap_type ca \
    --cmap_thresh ${cmap_thresh} \
    --train_tfrecord_fn ${data_dir}TFRecords/ \
    --valid_tfrecord_fn ${data_dir}TFRecords/ \
    --annot_fn ${annot_fn} \
    --test_list ${test_list} \
    --model_name ${model_name} \

    # --train_tfrecord_fn ${data_dir}${cmap_data}_${annot}_train \
    # --valid_tfrecord_fn ${data_dir}${cmap_data}_${annot}_valid \

    # --train_tfrecord_fn ${data_dir}${cmap_data}_train \
    # --valid_tfrecord_fn ${data_dir}${cmap_data}_valid \