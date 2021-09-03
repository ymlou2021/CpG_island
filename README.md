# CpG_island
compare CpG Islands found by python algorithms and UCSC

Task1 goal: 
According to Gardinaer-Garden and Frommer's paper, find the CpG islands in the first 1,000,000 base paris of human chromosome 20. Download the pG islands sequence from UCSC Genome Browser, and compare the CpG islands given by this two methods. Compare the overlap and create a confusion matrix to calculate the sensitivity , specigicity, false discovery rate and accuracy.
chromosome choosen: chr 20
human genome assembly: GRCh38

Workflow:
After running the python code, we get the cpg_py.bed file. Then the CpG islands (less than 1kbp apart) are merged: bedtools merge -i cpg_py_GRCh38.sort.bed -d 1000 > cpg_py_GRCh38.sort.merge.bed
Then the merge file is compared with CpG islands downloaded from UCSC: bedtools intersect -a cpg_ucsc.bed -b cpg_py_merge.bed -wo -f 0.4 | wc -l
