# Shands et al. (2023) _Phytophthora cinnamomi_ Genomics Materials
Shands et al., 2023 _(in prep)_

## Trimming Illumina Reads
The genomic Illumina reads were filtered using Trimmomatic v0.36. 
``` bash
exec java  -jar /opt/linux/centos/7.x/x86_64/pkgs/trimmomatic/0.36/trimmomatic-0.36.jar PE Pc2113_1.fq Pc2113_2.fq 2113_forward_paired.fq.gz 2113_forward_unpaired.fq.gz 2113_reverse_paired.fq.gz 2113_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
## *De novo* Genome Assmebly
Canu v1.7.1 (Koren et al. 2017), was used at default setting to generate both assemblies.
``` bash
canu \
-p 2113_Canu \
genomeSize=200m \
-pacbio-raw Pc2113_Combined_Subreads.fq
```

## Genome Assembly Correcting/Polishing
The resulting draft assemblies were corrected with PacBio reads iteratively three times using Racon v1.3.2 (Varser et al. 2017) follwed by polishing with the Illumina reads with Pilon v.1.22 (Walker et al. 2014) iteratively three times. 

Racon: 
``` bash

```


## Genome Size Estimation with Jellyfish & GenomeScope 
K-mer histograms were generated with trimmed Illumina reads using Jellyfish (v.2.3.0) (k-mer range 17-99 increments of 7). The respective histo outputs were uploaded to GenomeScope (http://qb.cshl.edu/genomescope/). 

Jellyfish count:
``` bash
for i in 17 25 33 41 49 57 65 73 81 89 97 99
do
  jellyfish count -m $i -s 1000000000  -t 4 -o Pc2113_k$i.jf 2113_forward_paired.fq 2113_reverse_paired.fq.gz
done
```
Jellyfish histo:
``` bash
for i in 17 25 33 41 49 57 65 73 81 89 97 99
do
  jellyfish histo -t 6 -o Pc2113_histo_k$i.txt Pc2113_k$i.jf
done
```

## Delimiting Genome into Gene-sparse & Gene-dense Regions 
This information can also be found in Supplementary Methods. To distinguish the gene-dense regions (GDR) and gene-sparse regions (GSR) of the genome we used methods described in Raffaele et al., 2010 and Rojas-Estevez et al., 2020. First, we simulated single-copy ‘core’ orthologs (n=2540) content determined by Orthofinder in GDRs and GSRs as percent of total genes belonging to each of these regions using values of the length ‘L’ of the FIRS between genes ranging from 100 bp to 5 Kb with 100bp increments. Genes with both FIRs greater than L were considered GSR genes and genes with both FIRs below L were considered GDR genes. Genes that had one FIR larger than L and the other lower than L were considered in-between, and genes with one FIR missing was considered not determined (ND). The core ortholog segregation rate was defined as the difference between the core ortholog content within the GDRs and GSRs, respectively. To determine the optimal L value that best fits the data, and that maximized the segregation rate and in which the percentage of core ortholog genes residing in GDR or in-between corresponded to at least 90%. 

Simulate_L.py was used to simulate L as follows for each isolate:
``` bash
python Simulate_L.py -i Pc2113_FIRs.csv -sco Pc2113_SCO.csv -s 100 -e 5100 -b 100
```
The results from this script were analzed in Excel to determine the optimal L-value. 


Next, we applied the optimal L-value to define the respective regions using Deliminate_Genome.py:
``` bash
python Deliminate_Genome.py -i Pc2113_SC_OG.csv -l 1100
```
