## build hominid phylogeny

# conda
conda activate muscle

# align seqs with mafft
mafft --auto data/hominid_cox1seqs.fasta > align_hominid_cox1seqs.fasta

# trim 
clipkit align_hominid_cox1seqs.fasta -o clipkit_hominid_cox1seqs.fasta

# build tree in iqtree
conda deactivate
micromamba activate iqtree  # ignore having conda and micromamba both running... rookie mistake
iqtree -s clipkit_hominid_cox1seqs.fasta --seqtype DNA -m TEST -bb 1000 --prefix hominid --seed 123

# open tree in figtree to root with Mangabey and export rooted tree 