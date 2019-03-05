rem ---- this script will reproduce the MPRAscore data from the paper.
rem ---- change -p switch to increase the p-value depth of the permutation testing (-p:500000 was used in the paper; but will take longer to compute).
rem ---- please note that the output file is not sorted by effect size ("score" column).

mprascore L363_barcode_map.txt L363_barcode_counts.txt -RNAcols:L363_RNA_1,L363_RNA_2,L363_RNA_3 -DNAcols:L363_DNA_1,L363_DNA_2,L363_DNA_3 L363_results.txt -p:1000
