#verilerin indirilmesi
bash setup/day2.sh

#dataların kontrolü
ls -lh data/filtered_gene_bc_matrices/hg19/
head -3 data/filtered_gene_bc_matrices/hg19/genes.tsv
head -3 data/filtered_gene_bc_matrices/hg19/barcodes.tsv
head -5 data/filtered_gene_bc_matrices/hg19/matrix.mtx

#github push
git add scripts/
git commit -m "Day 2: single-cell"
git push

