import pandas as pd

# Define terms to exclude
terms_to_remove = ['GO:0003674', 'GO:0005488', 'GO:0005515']  # molecular function, binding, protein binding

# Read usable_mf_terms.txt
with open('../APF/ALPHAFOLDpipeline/usable_mf_terms.txt', 'r') as f:
    mf_terms = sorted(set([line.strip() for line in f]))

# Read train_terms.tsv for discrepancy check
terms_df = pd.read_csv('../APF/Train/train_terms.tsv', sep='\t')
mf_terms_df = terms_df[(terms_df['aspect'] == 'MFO')]
pdb_to_mf_tsv = mf_terms_df.groupby('EntryID')['term'].apply(set).to_dict()

# Read usable_{train,val,test}_graphs.tsv to get splits, MF terms, and check discrepancies
splits = {'train': [], 'val': [], 'test': []}
pdb_to_mf = {}
discrepancies = []
num_ok = 0
for split in ['train', 'val', 'test']:
    df = pd.read_csv(f'../APF/ALPHAFOLDpipeline/usable_{split}_graphs.tsv', sep='\t')
    for _, row in df.iterrows():
        pdb_id = row['pdb_id']
        splits[split].append(pdb_id)
        # Get MF terms from label
        label = [int(x) for x in row['label'].split()]
        mf_terms_set = set(term for i, term in enumerate(mf_terms) if label[i] == 1 and term not in terms_to_remove)
        pdb_to_mf[pdb_id] = mf_terms_set
        # Reverse discrepancy check: label is ground truth
        tsv_terms = set(pdb_to_mf_tsv.get(pdb_id, set()))
        for i, term in enumerate(mf_terms):
            if (term in mf_terms_set and term not in tsv_terms) or (term not in mf_terms_set and term in tsv_terms):
                # discrepancies.append(f"Discrepancy in {split} for {pdb_id}: {term} (index={i}, label={term in mf_terms_set}, terms.tsv={term in tsv_terms})")
                discrepancies.append((f"Discrepancy in {split} for {pdb_id}: {term}", {"index": i, "label": term in mf_terms_set, "terms.tsv": term in tsv_terms}))
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
    if all([d['terms.tsv'] for s, d in discrepancies]):
        print("All discrepancies have a term in terms.tsv but not in the label. (Fine!)")
    elif all([not d['terms.tsv'] for s, d in discrepancies]):
        print("All discrepancies have a term in the label but not in terms.tsv.")
    else:
        print("Discrepancies have a term in both terms.tsv and the label.")
        for d in discrepancies:
            print(d[0], d[1])
    print(f"Number of correct labels: {num_ok}")
    print(f"Number of discrepancies: {len(discrepancies)}")
else:
    print("No discrepancies found.")