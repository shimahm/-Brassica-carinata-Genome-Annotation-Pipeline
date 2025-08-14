WIN=100000  # change your window here
> all_genes_around_snps.tsv
echo -e "snp_chr\tsnp_pos\tchr\tstart\tend\tgene_id\tlocus_tag\tortholog\tdescription\tGO\tGO_names" \
  > all_genes_around_snps.tsv

while read CHR POS; do
  python3 view_region_genes.py bcar_annotated.gff3 "$CHR" $((POS-WIN)) $((POS+WIN)) tmp.tsv
  awk -v c="$CHR" -v p="$POS" 'NR>1{print c"\t"p"\t"$0}' tmp.tsv >> all_genes_around_snps.tsv
done < snps.tsv
rm -f tmp.tsv

# Preview
head all_genes_around_snps.tsv

