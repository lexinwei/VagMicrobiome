beta=(bray-curtis jaccard weighted-unifrac unweighted-unifrac clr aitchison)
# t__
for m in ${beta[@]};do
  echo $m
  Rscript ~/workspace/vag/reference_code/MetaPhlAn-master/metaphlan/utils/calculate_diversity.R \
  -f ~/workspace/vag/metaphlan4_results/merged_abundance_table.txt \
  -o ~/workspace/vag/metaphlan4_results/diversity \
  -t /Users/xinwei/workspace/vag/reference_code/MetaPhlAn-master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk \
  -d beta -m $m -s t__ -p beta_T
done

alpha=(richness shannon simpson gini)
# t__
for m in ${alpha[@]};do
  echo $m
  Rscript ~/workspace/vag/reference_code/MetaPhlAn-master/metaphlan/utils/calculate_diversity.R \
  -f ~/workspace/vag/metaphlan4_results/merged_abundance_table.txt \
  -o ~/workspace/vag/metaphlan4_results/diversity \
  -t /Users/xinwei/workspace/vag/reference_code/MetaPhlAn-master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk \
  -d alpha -m $m -s t__ -p alpha_T
done