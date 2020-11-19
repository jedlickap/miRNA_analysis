# miRNA_analysis
miRNA_analysis

### 1. miRPara_outputParser

*This script filter the best candidates of miRNAs found by miRPara tool and produce \*.bed file and two fasta files with pre-miRNAs and mature miRNAs with correct headers like 'Alyr|scaffold_1|TE_132|copia_Ivana|noTSD'*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/miRPara_outputParser_v2.py \*_level_1.out \*_headerKeys.txt


### 2. miRPara_miRnaLtrRt_intersecs

*This script check the rough position of filtered miRNAs 'LTR-IN-LTR' and returns tables with all non-nested TEs with specific positions of miRNAs*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/miRPara_miRnaLtrRt_intersecs_v1.py *_TEs.fa *_level_1_filtered.bed *.gff

