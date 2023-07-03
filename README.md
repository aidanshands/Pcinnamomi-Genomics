# Shands et al. (2023) _Phytophthora cinnamomi_ Genomics Materials
_Shands et al., 2023 (in prep)_

## Trimming Illumina Reads
The genomic Illumina reads were filtered using Trimmomatic v0.36. 
``` bash
exec java  -jar /opt/linux/centos/7.x/x86_64/pkgs/trimmomatic/0.36/trimmomatic-0.36.jar PE Pc2113_1.fq Pc2113_2.fq 2113_forward_paired.fq.gz 2113_forward_unpaired.fq.gz 2113_reverse_paired.fq.gz 2113_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## Genome Size Estimation with Jellyfish & GenomeScope 
K-mer histograms were generated with trimmed Illumina reads using Jellyfish (v.2.3.0) (k-mer range 17-99 increments of 7) (Marçais and Kingsford, 2011.) The respective histo outputs were uploaded to GenomeScope (http://qb.cshl.edu/genomescope/). The same code was applied for Pc2109 with the respective genome assembly described in Shands _et al._ (2023).

**Jellyfish count:**
``` bash
for i in 17 25 33 41 49 57 65 73 81 89 97 99
do
  jellyfish count -m $i -s 1000000000  -t 4 -o Pc2113_k$i.jf 2113_forward_paired.fq 2113_reverse_paired.fq
  jellyfish histo -t 6 -o Pc2113_histo_k$i.txt Pc2113_k$i.jf
done
```

## *De novo* Genome Assmebly
Canu v1.7.1 (https://github.com/marbl/canu) (Koren _et al._ 2017), was used at default setting to generate both assemblies. The same code was applied for Pc2109 described in Shands _et al._ (2023).

``` bash
canu \
-p 2113_Canu \
genomeSize=200m \
-pacbio-raw Pc2113_Combined_Subreads.fq
```

## Genome Assembly Correcting/Polishing
The resulting draft assemblies were corrected with PacBio reads iteratively three times using Racon v1.3.2 (Varser _et al._ 2017) follwed by polishing with the Illumina reads with Pilon v.1.22 ([https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage](https://github.com/broadinstitute/pilon)) (Walker _et al._ 2014) iteratively three times. Minimap2 v2.10 (https://github.com/lh3/minimap2) (Li, 2018) was used to map the PacBio reads to the draft Canu assembly and the respective sam file was used for Racon. Bowtie2 v.2.3.5 (https://github.com/BenLangmead/bowtie2) (Langmead & Salzberg, 2012) and Samtools v. 1.17 (http://www.htslib.org/) (Li et al., 2009) were used for mapping and sam/bam processing prior to Pilon. The same code was applied for Pc2109 described in Shands _et al._ (2023).
**Racon Round 1:**
``` bash
minimap2 -ax map-pb -t8 Pc2113_Canu.contigs.fasta Pc2113_Combined_Subreads.fq > 2113_1.sam
~/racon/build/bin/racon -t 24 Pc2113_Combined_Subreads.fq 2113_1.sam Pc2113_Canu.contigs.fasta > Pc2113_racon1.fasta
```
**Racon Round 2:**
``` bash
minimap2 -ax map-pb -t8 Pc2113_racon1.fasta Pc2113_Combined_Subreads.fq > 2113_2.sam
~/racon/build/bin/racon -t 24 Pc2113_Combined_Subreads.fq 2113_2.sam Pc2113_racon1.fasta > Pc2113_racon2.fasta
```
**Racon Round 3:**
``` bash
minimap2 -ax map-pb -t8 Pc2113_racon2.fasta Pc2113_Combined_Subreads.fq > 2113_3.sam
~/racon/build/bin/racon -t 24 Pc2113_Combined_Subreads.fq 2113_3.sam Pc2113_racon2.fasta > Pc2113_racon3.fasta
```
**Pilon Round 1:**
``` bash
bowtie2-build Pc2113_racon3.fasta Pc2113_racon3
bowtie2 -x Pc2113_racon3 -1 2113_3_forward_paired.fq -2 2113_3_reverse_paired.fq | samtools view -Sb - | samtools sort -o Pc2113_1.bam -
samtools index -b Pc2113_1.bam
pilon --genome Pc2113_racon3.fasta --bam Pc2113_1.bam
# result: Pc2113_pilon1.fasta
```
**Pilon Round 2:**
``` bash
bowtie2-build Pc2113_pilon1.fasta Pc2113_pilon1
bowtie2 -x Pc2113_pilon1 -1 2113_3_forward_paired.fq -2 2113_3_reverse_paired.fq | samtools view -Sb - | samtools sort -o Pc2113_2.bam -
samtools index -b Pc2113_2.bam
pilon --genome Pc2113_pilon1.fasta --bam Pc2113_2.bam
# result: Pc2113_pilon2.fasta
```
**Pilon Round 3:**
``` bash
bowtie2-build Pc2113_pilon2.fasta Pc2113_pilon2
bowtie2 -x Pc2113_pilon2 -1 2113_3_forward_paired.fq -2 2113_3_reverse_paired.fq | samtools view -Sb - | samtools sort -o Pc2113_3.bam -
samtools index -b Pc2113_3.bam
pilon --genome Pc2113_pilon2.fasta --bam Pc2113_3.bam
# result: Pc2113_Polished.fasta
```

## Genome Assembly Purging
The resulted polished assemblies were subjected to haplotype purging to obtain a haploid genome assembly size that was consistent with our FCM estimations using Purge Haplotigs (https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) (Roach _et al._ 2018). QUAST v.4.6.3 (https://github.com/ablab/quast) (Gurevich _et al._ 2013) was used to evaluagte the genome assemblies. The same code was applied for Pc2109 described in Shands _et al._ (2023).
``` bash
# Minimap2
minimap2 -ax map-pb Pc2113_Polished.fasta Pc2113_Combined_Subreads.fq | samtools view -hF 256 - | samtools sort -@ 8 -o aligned.bam -T tmp
# Step 1
purge_haplotigs  readhist  -b aligned.bam  -g Pc2113_Polished.fasta  -t 8
# Step 2 the values here were used for Pc2109: -l 20  -m 83  -h 185
purge_haplotigs  contigcov  -i aligned.bam.gencov  -l 22  -m 97  -h 155  -o coverage_stats.csv
# Step 3
purge_haplotigs  purge  -g Pc2113_Polished  -c coverage_stats.csv  -o Pc2113_Purged  -a 79  -m 275  -t 16
# Quast
quast.py Pc2113_Purged.fasta
```

## Assessing Completeness via BUSCO
Benchmark Universal Single Copy Orthologs (BUSCO) v.3.0.2 (https://busco.ezlab.org/) (Simão _et al._ 2015) using the alveolata/stramenopiles database (orthoDB v.10) (Kriventseva _et al._ 2019) was used to assess completeness of the assemblies. The same code was applied for Pc2109 described in Shands _et al._ (2023).
``` bash
run_BUSCO.py -i Pc2113T1_genome.fasta \
-o Pc2113T1_AS \
-l alveolata_stramenophiles_ensembl \
-m geno \
-c 16
```

## Assessing Quality & Heterozygosity via K-mer Analysis Toolkit (KAT)
The quality and heterozygosity of the assemblies were assessed with KAT v2.3.4 (https://github.com/TGAC/KAT) (Mapleson _et al._ 2017). The same code was applied for Pc2109 described in Shands _et al._ (2023).

**KAT Comp**
``` bash
# Comp
kat comp -t 32 -o Pc2113 \
2113_forward_paired.fq 2113_reverse_paired.fq \
Pc2113T1_genome.fasta
```
**KAT Density**
``` bash
# Density 
kat comp -t 32 -n -o Pc2113 \
2113_forward_paired.fq 2113_reverse_paired.fq \
Pc2113T1_genome.fasta
```

## Nucmer/DNAdiff/Dotplot
Whole genome alignments of Pc2113 and Pc2109 was performed using Nucmer and differences were assessed using DNAdiff (https://mummer.sourceforge.net/) (Marçais _et al._ 2018). Whole-genome alignments for the construction of dot plots were generated using Minimap2 v2.10 (Li, 2018) and visualized using D-GENIES v1.3.1 (https://dgenies.toulouse.inra.fr/) (Cabanettes and Klopp, 2018). 
``` bash
# Nucmer
nucmer -t 32 -l 100 -c 500 --maxmatch Pc2113T1_genome.fasta Pc2109T1_genome.fasta -p 2113_vs_2109
# DNAdiff
dnadiff -p 2113_vs_2109 -d 2113_vs_2109.delta
# Minimap
minimap2 -t 32 -cx asm10 --cs Pc2113T1_genome.fasta Pc2109T1_genome.fasta > 2113_vs_2109.paf
```

## Ploidy
Ploidy analysis was performed on the _P. cinnamomi_ isolates and _P. infestans_ isolate 1306-C (Pan _et al._, 2108). Two methods were employed to determine ploidy; nQuire (Weiß _et al._, 2018) was used to estimate ploidy and the R (R Core Team 2020) package _vcfR_ (Knaus and Grünwald, 2017) was used to infer ploidy. For allele balance histograms generated in _vcfR_ I followed the tutorial by Brian J. Knaus and Niklaus J. Grünwald (https://knausb.github.io/vcfR_documentation/determining_ploidy_1.html). See CBS14422_Allele_Balance.R, Pc2113_Allele_Balance.R, Pc2109_Allele_Balance.R & Pi1306C_Allele_Balance.R. To run each of these methods, reads were filtered with fastq-mcf v1.05 (Aronesty, 2011) then BWA-mem (Li 2013) was used to map the filtered reads to their corresponding genome assemblies, and the respective BAM files were processed using Samtools (Li _et al._, 2009). Single-nucleotide polymorphisms (SNPs) were called from the processed BAM files using freebayes (Garrison and Marth, 2012). The same code was applied for all samples with the respective genome assembly described in Shands _et al._ (2023).

**Fastq-mcf**
``` bash
fastq-mcf \
-q 30 -D 150 -o Pc2113_1_filtered.fq -o Pc2113_2_filtered.fq \
30x_Illumina_Adapters.fasta \
Pc2113_1.fq \
Pc2113_2.fq
```

**BWA-mem**
``` bash
# performed the same for each isolate
bwa index -p Pc2113T1_genome.fasta
bwa mem -t 32 Pc2113T1_genome \
Pc2113_1_filtered.fq \
Pc2113_2_filtered.fq \
-R $(echo "@RG\tID:"${Pc2113}"\tSM:"${Pc2113}"\tLB:"${Pc2113}"\tPL:ILLUMINA")
| samtools view -bS - > Pc2113.bam
```

**Freebayes**
``` bash
freebayes -f Pc2113T1_genome.fasta \
--ploidy 2 \
CBS14422.bam \
Pc2109.bam \
Pc2113.bam | bgzip > Pc_Isolates_UNFILTERED.vcf.gz
```

**VCFtools**
``` bash
vcftools --gzvcf Pc_Isolates_UNFILTERED.vcf.gz \
--minQ 30 \
--thin 10 \
--out Pc_Isolates.F \
--recode-INFO-all \
--recode
```

**nQuire:**

The -c value for nQuire Create for each sample: Pc2113 (-c 40), CBS144.22 (-c 14), Pi1306-C (-c 13).
``` bash
# nQuire Create
~/nQuire/nQuire create \
-b Pc2109.bam \
-o Pc2109 \
-q 30 -c 31

# nQuire lrdmodel
~/nQuire/nQuire lrdmodel -t 8 \
Pc2109.bin > Pc2109_lrdmodel.txt

# nQuire Denoise
~/nQuire/nQuire denoise \
Pc2109.bin -o Pc2109_Denoised.bin

# nQuire lrdmodel on denoised bin
~/nQuire/nQuire lrdmodel -t 16 \
Pc2109_Denoised.bin > Pc2109_lrdmodel_Denoised.txt
```

## Orthofinder
Orthologs were identified using Orthofinder v.2.5.2 (Emms & Kelly, 2019).

``` bash
orthofinder -f Input_Fastas/
```

## Delimiting Genome into Gene-sparse & Gene-dense Regions 
This information can also be found in Supplementary Methods. To distinguish the gene-dense regions (GDR) and gene-sparse regions (GSR) of the genome we used methods described in Raffaele et al., 2010 and Rojas-Estevez et al., 2020. First, we simulated single-copy ‘core’ orthologs (n=2540) content determined by Orthofinder in GDRs and GSRs as percent of total genes belonging to each of these regions using values of the length ‘L’ of the FIRS between genes ranging from 100 bp to 5 Kb with 100bp increments. Genes with both FIRs greater than L were considered GSR genes and genes with both FIRs below L were considered GDR genes. Genes that had one FIR larger than L and the other lower than L were considered in-between, and genes with one FIR missing was considered not determined (ND). The core ortholog segregation rate was defined as the difference between the core ortholog content within the GDRs and GSRs, respectively. To determine the optimal L value that best fits the data, and that maximized the segregation rate and in which the percentage of core ortholog genes residing in GDR or in-between corresponded to at least 90%. The same code was applied for Pc2109 described in Shands _et al._ (2023).

Simulate_L.py was used to simulate L as follows for each isolate:
``` bash
python Simulate_L.py -i Pc2113_FIRs.csv -sco Pc2113_SCO.csv -s 100 -e 5100 -b 100
```
The results from this script were analzed in Excel to determine the optimal L-value. 

Next, we applied the optimal L-value to define the respective regions using Deliminate_Genome.py:
``` bash
python Deliminate_Genome.py -i Pc2113_SC_OG.csv -l 1100
```

## Secretome & RXLR Effectors
Proteome for each _P. cinnamomi_ isolate was scanned for signal peptide (SP) presence using 
SignalP (v.5.0) (Armenteros _et al._ 2019). The resulting proteins containing SP sequence were subjected to TMHMM (v.2.0) analyses (Möller _et al._, 2001) to identify transmembrane domains (TMDs). TMD-containing proteins were removed from the secretome protein dataset. Secretomes were subjected to EffectorP v3.0 (Sperschneider & Dodds, 2021) to predict apoplastic, cytoplasmic, and dual-localized effectors. RXLR effectors were predicted from the secretome using regular expression (REGEX) searches (FindRXLRs.py). The same code was applied for Pc2109 described in Shands _et al._ (2023).

**Signal P**
``` bash
signalp -fasta Pc2113T1_genome.pep.fasta -org euk -format short -prefix Pc2113_Short
```

**TMHMM**
``` bash
tmhmm --short Pc2113T1_SP5.fasta > Pc2113T1_SP5_TMHMM.txt
```

**Effector P 3.0**
``` bash
python EffectorP.py -i Pc2113T1_SP5.fasta > Pc2113_EffectorP_out.txt
```

**Finding RXLRs**
``` bash
python FindRXLRs.py -i Pc2113T1_SP5_noTMHMM.fasta
```

**IQ-Tree2**
``` bash
# Muscle v.3.8.425 alignment 
muscle -in All_RXLRs.fasta -out RXLR_Alignment.fasta

# IQ-Tree v. 2.1.3 Modelfinder
iqtree2 -s RXLR_Alignment.fasta -mem 32G -T AUTO -m MFP

# IQ-Tree v. 2.1.3 
iqtree2 -s RXLR_Alignment.fasta -mem 32G -T AUTO -m VT+F+R7 -bb 1000 -nm 5000

```

