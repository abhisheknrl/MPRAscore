# MPRAscore
Robust and non-parametric analysis of massively parallel reporter assays (MPRA)

MPRAscore infers allele-specific effects on transcription from MPRA data. MPRAscore uses a weighted, variance-regularized method to calculate variant effect sizes robustly, and a permutation approach to test for significance without assuming normality or independence.


Running MPRAscore:
MPRAscore is implemented in C++. It can be compiled from the source code. Precompiled executables for Windows, Linux, and Mac are included in binaries/.

For help:
mprascore -?


Example:
Example files are included in Example/. To run MPRAscore on example data:

mprascore L363_barcode_map.txt L363_barcode_counts.txt -RNAcols:L363_RNA_1,L363_RNA_2,L363_RNA_3 -DNAcols:L363_DNA_1,L363_DNA_2,L363_DNA_3 L363_results.txt -p:1000


