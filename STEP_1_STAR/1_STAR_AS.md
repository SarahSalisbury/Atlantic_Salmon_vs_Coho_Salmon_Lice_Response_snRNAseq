# 1_STAR_AS

This document outlines how to use STAR to index the Atlantic Salmon Genome (V3) from ENSEMBL and then map your single cell RNAseq libraries to this indexed genome.

## Get Genome Files
Download the ENSEMBL version of the V3 genome.

Download the fasta file (Salmo_salar.Ssal_v3.1.dna_rm.toplevel.fa.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/.

Download the gff file (Salmo_salar.Ssal_v3.1.106.gff3.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-106/gff3/salmo_salar/.

Gunzip the files
```
gunzip *.gz
```

## Get Mitochondrial DNA Information
Ok but before we can run any indexing, note that there is no mtDNA in the V3 version!

So what we need is to add the mitochondria dna to our fasta and their annotations to our gff3 file.

We're going to use the V2 ENSEMBL version of the genome to do this.

So, download the fasta and gff of the ENSEMBL version V2.

Download the fasta file (Salmo_salar.ICSASG_v2.dna_rm.toplevel.fa.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-105/fasta/salmo_salar/dna/.

Download the gff file (Salmo_salar.ICSASG_v2.105.gff3.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-105/gff3/salmo_salar/.

Gunzip the files
```
gunzip *.gz
```

### Append mitochondrial genome to nuclear genome for .gff file

So for the V2 ENSEMBL Atlantic Salmon genome, the chromosome for the mitochondria is "MT".

So we can just find all the lines with MT (mtDNA code) at the beginning of the line in the V2 genome .gff, note that there is no header when you do this (so you don't have to worry about removing it before appending it to the .gff of the nuclear V3 genome.
```
grep '^MT*' Salmo_salar.ICSASG_v2.105.gff3 > ./mtDNAV2nohead.gff
```

Maybe you'd like to know what the genes are called in the mitochondria (we will need these later in Seurat when we want to screen out mitochondrial genes as a quality control measure).
So what we'll do first is grab the lines which have a gene feature (have "protein_coding" within them), next of these grep the lines that have "ID=gene", then we'll use sed to save everything in the line which is after (any number of characters - .(asterisk))"gene=", next we'll use sed to replace anything following a ";" with nothing (i.e. deleting it), then we'll sort the remaining values alphabetically, and then just get the uniq values
```
grep "protein_coding" mtDNAV2nohead.gff | grep "ID=gene" |sed 's/.*Name=\(.*\)/\1/' | sed 's/;.*//'| sort | uniq


ATP6
ATP8
COX1
COX2
COX3
CYTB
ND1
ND2
ND3
ND4
ND4L
ND5
ND6

```

Maybe you also want to know the gene IDs. This is similar to what we did above, just get rid of everything before ID=gene: (inclusive).
```
grep "protein_coding" mtDNAV2nohead.gff | grep "ID=gene" |sed 's/.*ID=gene:\(.*\)/\1/' | sed 's/;.*//'| sort | uniq

ENSSSAG00000000007
ENSSSAG00000000011
ENSSSAG00000000017
ENSSSAG00000000020
ENSSSAG00000000022
ENSSSAG00000000023
ENSSSAG00000000024
ENSSSAG00000000026
ENSSSAG00000000028
ENSSSAG00000000029
ENSSSAG00000000033
ENSSSAG00000000034
ENSSSAG00000000036
```

Now make a copy of the V3 nuclear genome gff
```
cp Salmo_salar.Ssal_v3.1.106.gff3 Salmo_salar.Ssal_v3.1.106_MT.gff
```

Note that before doing the next step you might have to change permissions of GCF_905237065.1_Ssal_v3.1_genomic_MT.gff to make it writable ```+w```
```
chmod +w Salmo_salar.Ssal_v3.1.106_MT.gff
```

Append mtdna to gff V3 nuclear genome
```
cat mtDNAV2nohead.gff >> Salmo_salar.Ssal_v3.1.106_MT.gff
```

### Append mitochondrial genome to nuclear genome for .fna file

Isolate the mtDNA fasta read from the v2 genome:
```
sed -n -e '/^>MT/,/>/ p' Salmo_salar.ICSASG_v2.dna_rm.toplevel.fa | sed -e '$d' > MTDNAV2.fna
```
NOTES:
- ```-e```, needed because you're running sed twice, you're providing it with multiple commands so you need to put the ```-e``` to specify each time you run sed.
- ```$d```, deletes the last line, specifically: ```$``` indicates last line, ```d``` indicates delete

Ok now make a copy of the V3 nuclear genome .fna file.
```
cp Salmo_salar.Ssal_v3.1.dna_rm.toplevel.fa Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna
```

Note that before doing the next step you might have to change permissions of Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna to make it writable ```+w```
```
chmod +w Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna
```

Append mtdna .fna file to fasta for the rest of the V3 nuclear genome. Note that there's no headers or anything in a .fna file, so no need to worry about appending one file to another.
```
cat MTDNAV2.fna >> Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna
```

You now have a .fna and .gff file containing both the nuclear V3 genome and the V2 mitochondrial genome!

## Prepare to Index!

### 1. Convert GFF to GTF file
Transform the gff genome annotation file to gft format for STARsolo using gffread:
```
gffread Salmo_salar.Ssal_v3.1.106_MT.gff -T -o Salmo_salar.Ssal_v3.1.106_MT.gtf
```
NOTES:
- ```-T```, specifies that you want it in GTF2 version, rather than default GFF3 version
- ```-o```, specifies you want to name the output file something specific (as specified after the flag).

### 2. Get Whitelist
We need the whitelist prior to indexing, this is the list of cell barcodes and is specific to the chemistry of the 10X kit you  used. https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-

We used the v3 kits so we need the file: ```3M-february-2018.txt```. Which can be found here: https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

### 3. Make output folder for indexing step
So we need to make a folder that we will then output our STAR outputs, this is specified in the indexing code below by ```--genomeDir star_genome_index```
Make star_genome_index
```
mkdir star_genome_index
```

## Indexing
```
STAR --limitGenomeGenerateRAM 586158178570 --runThreadN 10 --runMode genomeGenerate --genomeDir star_genome_index --genomeFastaFiles Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna --sjdbGTFfile Salmo_salar.Ssal_v3.1.106_MT.gtf
```
NOTES:

```--limitGenomeGenerateRAM``` indicates total amount of RAM STAR can use

```--runThreadN``` option tells you how many threads to have

```--runMode``` genomeGenerate option directs STAR to run genome indices generation job

```--genomeDir``` specifies path to the directory (henceforth called "genome directory” where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to have writing permissions. The file system needs to have at least 100GB of disk space available for a typical mammalian genome. It is recommended to remove all files from the genome directory before running the genome generation step. This directory path will have to be supplied at the mapping step to identify the reference genome. SO MAKE SURE YOU HAVE MADE THE SPECIFIED GENOME DIRECTORY BEFORE RUNNING

```--genomeFastaFiles``` specifies your fasta file

```--sjdbGTFfile``` /path/to/annotations.gtf

## Mapping
For each of our libraries we ran a script like the following:
```
STAR --genomeDir refv3.2/star_genome_index --readFilesIn SampleX_R2.fastq.gz SampleX_R1.fastq.gz --soloMultiMappers Unique EM --soloBarcodeReadLength 28 --soloType CB_UMI_Simple --soloUMIlen 12 --soloCBwhitelist refv3.2/3M-february-2018.txt --soloFeatures GeneFull --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 4 --outMultimapperOrder Random --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10  --readFilesCommand zcat --outSAMtype BAM Unsorted
```
NOTES:
```--genomeDir``` is the location where your indexed genome is

```--readFilesIn``` is your raw library reads (read2 then read1)

```--soloMultiMappers Unique EM``` is the counting method for reads mapping to multiple genes. So "Unique" is the default, where you count only reads that map to unique genes. "EM" specifies that UMIs which map to multiple genes are "counted" distributed across all the genes they map to using the Expectation Maximization algorithm, so basically you're splitting up your UMI "count" over many possible genes (so your final count table might have non-integer values).

```--soloBarcodeReadLength``` is telling star how long you expect your read to be, the default is the length of your UMI + the length of your barcode (so 12 bp for UMI + 16bp for barcode =28 bp)

```--soloType CB_UMI_Simple``` indicates that you have 10X Chromium libraries (so star should look for one UMI and one cell barcode)

```--soloUMIlen``` the length of the UMI, the default here is 10 but the new chemistry of 10X uses 12bp UMI, so this must be specified!

```--soloCBwhitelist``` this is our list of cell barcodes (CB) we expect (this is the file we got from Chromium earlier)

```--soloFeatures GeneFull``` so this is specifying which genomic features we want associated with our UMI and cell barcodes. "GeneFull" indicates that we want full gene (pre-mRNA) info, including if reads are in genes' exonns and introns.

```--clipAdapterType CellRanger4``` so this just specifies that we want star to cut out our adapters (5p and 3p) as is also done in CellRanger4, Utilizes Opal package by Martin ˇSoˇsi ́c: https://github.com/Martinsos/opal.

```--outFilterScoreMin``` integer: alignment will be output only if its score is higher than or equal to this value, default = 0

```--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts``` multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches AND pseudocounts of 1 are added to all whitelist barcodes AND  multimatching to whitelist is allowed for CBs with N-bases. This option matches best with CellRanger >=3.0.0

```--soloUMIfiltering MultiGeneUMI_CR``` remove UMIs with N and homopolymers AND remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0

```--soloUMIdedup 1MM_CR``` CellRanger2-4 algorithm for 1MM UMI collapsing

```--runThreadN``` defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node

```--outMultimapperOrder Random``` order of multimapping alignments in the output files random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases.

```--outFilterMultimapScoreRange 1``` integer: the score range below the maximum score for multimapping alignments, default is 1

```--outFilterMultimapNmax 2``` int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as ”mapped to too many loci” in the Log.final.out. Default is 10

```--outFilterMismatchNmax 10``` integer: alignment will be output only if it has no more mismatches than this value. Default is 10

```--readFilesCommand zcat``` command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout, zcat will uncompress .gz files

```--outSAMtype BAM Unsorted``` output BAM without sorting

SANITY CHECKS:
IS YOUR ```--genomeDir``` THE SAME AS THE ONE YOU JUST MADE IN THE INDEX STEP?

IS ```--readFilesIn``` ASSOCIATED WITH THE CORRECT FILES AND IS READ2 BEFORE READ1?

> [!NOTE]  
> Even though mapping included ```--soloMultiMappers Unique EM``` we only used the matrix file generated from the "unique" mapping for downstream analysis in Seurat (i.e., ```matrix.mtx```) not the one output with unique and multimapped reads (```UniqueAndMult-EM.mtx```). Therefore the following parameters were irrelevant ```--outMultimapperOrder Random``` ```--outFilterMultimapScoreRange 1``` ```--outFilterMultimapNmax 10```. 
