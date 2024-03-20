# 2_STAR_CO

This document outlines how to use STAR to index the Coho Salmon Genome (V2) from ENSEMBL and then map your single cell RNAseq libraries to this indexed genome.

## Get Genome Files
Download the ENSEMBL version of the V3 genome.

Download the fasta file (Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel.fa.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-106/fasta/oncorhynchus_kisutch/dna/.

Download the gff3 file (Oncorhynchus_kisutch.Okis_V2.106.gff3.gz) and checksums (CHECKSUMS) from http://ftp.ensembl.org/pub/release-106/gff3/oncorhynchus_kisutch/.

Gunzip the files
```
gunzip *.gz
```

## Get Mitochondrial DNA Information
Ok but before we can run any indexing, we need the MTDNA!

Because there is no ensembl version of the mtDNA for Coho, we're going to use the mtDNA from NCBI associated with this version of the coho genome.

Specifically, we'll use the mitochondrial genome associated with the RefSeq version of the ENSEMBL annotation of the Okis V2 Genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_009263.1.

### Clean up the GFF3 file

Delete all lines that start with "#"
https://stackoverflow.com/questions/8206280/delete-all-lines-beginning-with-a-from-a-file
```
sed '/^#/d' sequence.gff3 > sequencenohead.gff3
```

### Append mitochondrial genome to nuclear genome for .gff file

So for the Okis_V2 genome, the chromosome for the mitochondria is "NC_009263.1".

Maybe you'd like to know what the genes are called in the mitochondria (we will need these later in Seurat when we want to screen out mitochondrial genes as a quality control measure).
So what we'll do first is grab the lines which have a gene feature (have "protein_coding" within them), next of these grep the lines that have "ID=gene", then we'll use sed to save everything in the line which is after (any number of characters - .(asterisk))"gene=", next we'll use sed to replace anything following a ";" with nothing (i.e. deleting it), then we'll sort the remaining values alphabetically, and then just get the uniq values
```
grep "protein_coding" sequencenohead.gff3 | grep "ID=gene" |sed 's/.*Name=\(.*\)/\1/' | sed 's/;.*//'| sort | uniq

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

Maybe you also want to know the gene IDs. So this is similar to what we did above, just get rid of everything before ID=gene: (inclusive).
```
grep "protein_coding" sequencenohead.gff3 | grep "ID=gene" |sed 's/.*ID=gene-\(.*\)/\1/' | sed 's/;.*//'| sort | uniq

KEG53_p01
KEG53_p02
KEG53_p03
KEG53_p04
KEG53_p05
KEG53_p06
KEG53_p07
KEG53_p08
KEG53_p09
KEG53_p10
KEG53_p11
KEG53_p12
KEG53_p13
```

Ok but there is something funny going on with the annotations for the mtDNA genes! Only CDS and gene is listed for each, there is no RNA or exon. The latter in particular is a problem since STAR maps reads to exons! Ok so what we're going to do is change the CDS rows to exon rows. We're going to try to be as minimally invasive as possible by first changing all instances of "CDS" to "exon" and "cds" to "exon".

You can check to make sure that "CDS" and "cds" aren't in any gene names or something that we wouldn't want to change with the following lines of code:
```
grep 'CDS' sequencenohead.gff3
grep 'cds' sequencenohead.gff3
```

Ok now that we're confident that we can change "CDS" and "cds" without messing up our file, let's substitute "CDS" and "cds" for "exon" using sed and a global (g) substitution.
```
sed 's:CDS:exon:g' sequencenohead.gff3 > sequencenoheadCDStoexon.gff3
sed 's:cds:exon:g' sequencenoheadCDStoexon.gff3 > sequencenoheadCDSandcdstoexon.gff3
```

We would also like to know the gene IDs of all the genes including the tRNAs and rRNAs in the mtDNA (as these may also be amplified in our libraries). Unlike the Atlantic Salmon mtDNA, here even the tRNAs and rRNAs have "Names", so we can just use the above code and look for all genes with "Names" (which we'll need when filtering out mtDNA features during the Seurat analysis).
```
grep "ID=gene" sequencenoheadCDSandcdstoexon.gff3 | sed 's/.*Name=\(.*\)/\1/' | sed 's/;.*//'| sort | uniq

ATP6
ATP8
COX1
COX2
COX3
CYTB
KEG53_r01
KEG53_r02
KEG53_t01
KEG53_t02
KEG53_t03
KEG53_t04
KEG53_t05
KEG53_t06
KEG53_t07
KEG53_t08
KEG53_t09
KEG53_t10
KEG53_t11
KEG53_t12
KEG53_t13
KEG53_t14
KEG53_t15
KEG53_t16
KEG53_t17
KEG53_t18
KEG53_t19
KEG53_t20
KEG53_t21
KEG53_t22
ND1
ND2
ND3
ND4
ND4L
ND5
ND6

```
Please note that we will need the gene symbols for the named mtDNA genes (e.g., ND1) AND the gene ids for the unnamed genes in our downstream analysis in Seurat. So keep this information handy!

Now make a copy of the V3 nuclear genome gff
```
cp Oncorhynchus_kisutch.Okis_V2.106.gff3 Oncorhynchus_kisutch.Okis_V2.106_MT.gff
```

Note that before doing the next step you might have to change permissions of GCF_905237065.1_Ssal_v3.1_genomic_MT.gff to make it writable ```+w```
```
chmod +w Oncorhynchus_kisutch.Okis_V2.106_MT.gff
```

Append mtdna to gff V3 nuclear genome
```
cat sequencenoheadCDSandcdstoexon.gff3 >> Oncorhynchus_kisutch.Okis_V2.106_MT.gff
```

### Append mitochondrial genome to nuclear genome for .fna file

Since the sequence.fasta file only contains the mtDNA, we don't need to filter it. Notice though that the length of the line in the sequence.fasta is longer than that in Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel.fa (it's 70 bp long in former compared to 60 bp in the latter). Note that we didn't have to deal with this for the Atlantic Salmon genome, because we used the mtDNA from the ENSEMBL V2 genome and the width of that fasta file was 60 bp as was the width of the ENSEMBL V3 genome!

So first save the first line of the sequence.fasta file:
```
head -n 1 sequence.fasta > sequencehead.fasta
```

Now delete first line of the sequence.fasta file:
```
tail -n +2 sequence.fasta > sequencenohead.fasta
```

Now remove all breaks in the fasta:
https://stackoverflow.com/questions/3134791/how-do-i-remove-newlines-from-a-text-file
```
tr -d '\n' < sequencenohead.fasta > sequencenoheadnofastanobreaks.fasta
```

Now introduce a break every 60 bp:
https://stackoverflow.com/questions/1187078/how-to-insert-a-new-line-character-after-a-fixed-number-of-characters-in-a-file
```
sed -e "s/.\{60\}/&\n/g" < sequencenoheadnofastanobreaks.fasta > sequencenoheadnofastanobreaks60bpsplit.fa
```

Now compend the new sequence to the original heading:
```
cat sequencenoheadnofastanobreaks60bpsplit.fa >> sequencehead.fasta
```

Ok now make a copy of the nuclear genome .fna file.
```
cp Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel.fa  Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel_MT.fa 
```

Note that before doing the next step you might have to change permissions of Salmo_salar.Ssal_v3.1.dna_rm.toplevel_MT.fna to make it writable ```+w```
```
chmod +w Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel_MT.fa
```

Append mtdna .fna file to fasta for the rest of the V3 nuclear genome. Note that there's no headers or anything in a .fna file, so no need to worry about appending one file to another.
```
cat sequencehead.fasta >> Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel_MT.fa
```

You now have a .fa and .gff file containing both the nuclear genome and the mitochondrial genome!

## Prepare to Index!

### 1. Convert GFF to GTF file
Transform the gff genome annotation file to gft format for STARsolo using gffread:
```
gffread Oncorhynchus_kisutch.Okis_V2.106_MT.gff -T -o Oncorhynchus_kisutch.Okis_V2.106_MT.gtf
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
STAR --limitGenomeGenerateRAM 1200000000000 --runThreadN 10 --runMode genomeGenerate --genomeDir star_genome_index_moremem --genomeFastaFiles Oncorhynchus_kisutch.Okis_V2.dna_rm.toplevel_MT.fa --sjdbGTFfile Oncorhynchus_kisutch.Okis_V2.106_MT.gtf
```
NOTES:

```--limitGenomeGenerateRAM``` option indicates total amount of RAM STAR can use

```--runThreadN``` option tells you how many threads to have

```--runMode``` genomeGenerate option directs STAR to run genome indices generation job

```--genomeDir``` specifies path to the directory (henceforth called "genome directory” where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to have writing permissions. The file system needs to have at least 100GB of disk space available for a typical mammalian genome. It is recommended to remove all files from the genome directory before running the genome generation step. This directory path will have to be supplied at the mapping step to identify the reference genome. SO MAKE SURE YOU HAVE MADE THE SPECIFIED GENOME DIRECTORY BEFORE RUNNING

```--genomeFastaFiles``` specifies your fasta file

```--sjdbGTFfile``` /path/to/annotations.gtf

## Mapping
For each of our libraries we ran a script like the following:
```
STAR --genomeDir ../../refcoho/star_genome_index_moremem --readFilesIn ../../upto33/16-fhf_R2_001.fastq.gz ../../upto33/16-fhf_R1_001.fastq.gz --soloMultiMappers Unique EM --soloBarcodeReadLength 28 --soloType CB_UMI_Simple --soloUMIlen 12 --soloCBwhitelist ../../refv3.2/3M-february-2018.txt --soloFeatures GeneFull --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 4 --outMultimapperOrder Random --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10  --readFilesCommand zcat --outSAMtype BAM Unsorted
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

```--outFilterMultimapNmax 10``` int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as ”mapped to too many loci” in the Log.final.out. Default is 10

```--outFilterMismatchNmax 10``` integer: alignment will be output only if it has no more mismatches than this value. Default is 10

```--readFilesCommand zcat``` command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout, zcat will uncompress .gz files

```--outSAMtype BAM Unsorted``` output BAM without sorting

SANITY CHECKS:
IS YOUR ```--genomeDir``` THE SAME AS THE ONE YOU JUST MADE IN THE INDEX STEP?

IS ```--readFilesIn``` ASSOCIATED WITH THE CORRECT FILES AND IS READ2 BEFORE READ1?

> [!NOTE]
> Even though mapping included ```--soloMultiMappers Unique EM``` we only used the matrix file generated from the "unique" mapping for downstream analysis in Seurat (i.e., ```matrix.mtx```) not the one output with unique and multimapped reads (```UniqueAndMult-EM.mtx```). Therefore the following parameters were irrelevant ```--outMultimapperOrder Random``` ```--outFilterMultimapScoreRange 1``` ```--outFilterMultimapNmax 10```. 
