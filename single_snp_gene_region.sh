# set your SNP and window
CHR="CM081020.1"; POS=31630794; WIN=100000   # ±100 kb example
python3 view_region_genes.py \
  bcar_annotated.gff3 "$CHR" $((POS-WIN)) $((POS+WIN)) \
  region_all.tsv

# Preview
head region_all.tsv

