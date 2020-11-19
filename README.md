# miRNA_analysis
miRNA_analysis

### 1. miRPara_outputParser

*This script filters the best candidates of miRNAs found by miRPara tool and produce \*.bed file and two fasta files with pre-miRNAs and mature miRNAs with correct headers like 'Alyr|scaffold_1|TE_132|copia_Ivana|noTSD'*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/miRPara_outputParser_v2.py \*_level_1.out \*_headerKeys.txt


### 2. miRPara_miRnaLtrRt_intersecs

*This script checks the rough position of filtered miRNAs 'LTR-IN-LTR' and returns tables with all non-nested TEs with specific positions of miRNAs*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/miRPara_miRnaLtrRt_intersecs_v1.py *_TEs.fa *_level_1_filtered.bed *.gff


### 3. nesterMergedGffOUTsParser_allProtDom

*This script filters only thouse TEs with all protein domains detected and all therespective miRNAs within these TEs. The gff file with those TEs and three fasta files with those TEs, miRNAs and pre-miRNAs within them plus miRNA bed file are returned*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/nesterMergedGffOUTsParser_allProtDom.py *_miRNA.bed *.gff *_TEs.fa *_miRna.fa *_preMiRna.fa


### 4. exactMiRnaPosInTe

*This script scritp checks the accurate position of filtered miRNAs and then returns output tables '\*_miRNA_accuratePosInTe.tsv'*

> $ python3 /home/pavel/Documents/Soft/miRNA_analysis/exactMiRnaPosInTe.py *_filtered_allProtDomTes.bed *_filtered_allProtDomTes.gff
