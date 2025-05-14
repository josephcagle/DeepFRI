import os
import numpy as np
import joblib
from Bio import SeqIO
from torch_geometric.data import Data
from tqdm import tqdm

# Paths
graph_dir = "../APF/AlphaFold Data/graphs/"
fasta_file = "../APF/Train/train_sequences.fasta"
output_dir = "preprocessing/data2/annot_pdb_chains_npz/"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load fasta sequences
fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Counter for mismatches
mismatch_count = 0

# Process each joblib file
files = [x for x in os.listdir(graph_dir) if x.endswith("_0.graph.joblib")]
for filename in tqdm(files, total=len(files), smoothing=0):
    uniprot_id = filename.split("_")[0]
    
    # Load graph
    graph = joblib.load(os.path.join(graph_dir, filename))
    
    # Get sequence
    seqres = fasta_dict.get(uniprot_id)
    if seqres is None:
        raise ValueError(f"No sequence found for {uniprot_id}")
    
    # Check dimensions
    n = graph.num_nodes
    if n != len(seqres):
        mismatch_count += 1
        continue
        
    # Create adjacency matrix
    c_alpha = np.full((n, n), 13, dtype=np.float32)
    edge_index = graph.edge_index.numpy()
    c_alpha[edge_index[0], edge_index[1]] = 1
    c_alpha[edge_index[1], edge_index[0]] = 1  # Ensure symmetry
    
    # Save to npz
    output_path = os.path.join(output_dir, f"{uniprot_id}-A.npz")
    np.savez(output_path, C_alpha=c_alpha, C_beta=[], seqres=seqres)

print(f"Number of mismatches: {mismatch_count}")
