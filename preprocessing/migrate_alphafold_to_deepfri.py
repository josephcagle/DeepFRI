import pandas as pd

# Define root terms to exclude
terms_to_remove = ['GO:0003674', 'GO:0005488', 'GO:0005515'] # molecular function, binding and protein binding
# terms_to_remove = set()

# Read usable_mf_terms.txt
with open('../APF/ALPHAFOLDpipeline/usable_mf_terms.txt', 'r') as f:
    mf_terms = sorted(set([line.strip() for line in f]))

# Read train_terms.tsv and filter for MFO
terms_df = pd.read_csv('../APF/Train/train_terms.tsv', sep='\t')
mf_terms_df = terms_df[(terms_df['aspect'] == 'MFO')]
pdb_to_mf = mf_terms_df.groupby('EntryID')['term'].apply(set).to_dict()

# Read usable_{train,val,test}_graphs.tsv to get splits and check discrepancies
splits = {'train': [], 'val': [], 'test': []}
discrepancies = []
num_ok = 0
for split in ['train', 'val', 'test']:
    df = pd.read_csv(f'../APF/ALPHAFOLDpipeline/usable_{split}_graphs.tsv', sep='\t')
    for _, row in df.iterrows():
        pdb_id = row['pdb_id']
        splits[split].append(pdb_id)
        # Check for discrepancies
        label = [int(x) for x in row['label'].split()]
        expected_terms = set(pdb_to_mf.get(pdb_id, set()))
        for i, term in enumerate(mf_terms):
            if (term in expected_terms and label[i] == 0) or (term not in expected_terms and label[i] == 1):
                discrepancies.append(f"Discrepancy in {split} for {pdb_id}: {term} (index={i}, label={label[i]}, terms.tsv={term in expected_terms})")
            else:
                num_ok += 1

# Write split files
for split in ['train', 'val', 'test']:
    with open(f'./preprocessing/data2/alphafold_{split}.txt', 'w') as f:
        for pdb_id in splits[split]:
            f.write(f"{pdb_id}-A\n")

# Write alphafold_annot.tsv
with open('./preprocessing/data2/alphafold_annot.tsv', 'w') as f:
    f.write("### GO-terms (molecular_function)\n")
    f.write('\t'.join([term for term in mf_terms if term not in terms_to_remove]) + '\n')
    f.write("### GO-names (molecular_function)\n\n")
    f.write("### GO-terms (biological_process)\n\n")
    f.write("### GO-names (biological_process)\n\n")
    f.write("### GO-terms (cellular_component)\n\n")
    f.write("### GO-names (cellular_component)\n\n")
    f.write("### PDB-chain\tGO-terms (molecular_function)\tGO-terms (biological_process)\tGO-terms (cellular_component)\n")
    for pdb_id in sorted(set(splits['train'] + splits['val'] + splits['test'])):
        mf_terms_set = pdb_to_mf.get(pdb_id, set())
        mf_terms_str = ','.join(sorted([term for term in mf_terms_set if term not in terms_to_remove])) if mf_terms_set else ''
        f.write(f"{pdb_id}-A\t{mf_terms_str}\t\t\n")

# Print discrepancies
if discrepancies:
    print("Discrepancies found:")
    for d in discrepancies:
        print(d)
    print(f"Number of correct labels: {num_ok}")
    print(f"Number of discrepancies: {len(discrepancies)}")
else:
    print("No discrepancies found.")