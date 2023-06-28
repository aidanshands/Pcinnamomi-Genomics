# Shands et al. (2023) _Phytophthora cinnamomi_ Genomics Materials
Shands et al., 2023 _(in prep)_

## Delimiting Genome into Gene-sparse & Gene-dense Regions 
This information can also be found in Supplementary Methods. To distinguish the gene-dense regions (GDR) and gene-sparse regions (GSR) of the genome we used methods described in Raffaele et al., 2010 and Rojas-Estevez et al., 2020. First, we simulated single-copy ‘core’ orthologs (n=2540) content determined by Orthofinder in GDRs and GSRs as percent of total genes belonging to each of these regions using values of the length ‘L’ of the FIRS between genes ranging from 100 bp to 5 Kb with 100bp increments. Genes with both FIRs greater than L were considered GSR genes and genes with both FIRs below L were considered GDR genes. Genes that had one FIR larger than L and the other lower than L were considered in-between, and genes with one FIR missing was considered not determined (ND). The core ortholog segregation rate was defined as the difference between the core ortholog content within the GDRs and GSRs, respectively. To determine the optimal L value that best fits the data, and that maximized the segregation rate and in which the percentage of core ortholog genes residing in GDR or in-between corresponded to at least 90%. 

Simulate_L.py was used to simulate L as follows for each isolate:
``` bash
python Simulate_L.py -i Pc2109_FIRs.csv -sco Pc2109_SCO.csv -s 100 -e 5100 -b 100

python Simulate_L.py -i Pc2113_FIRs.csv -sco Pc2113_SCO.csv -s 100 -e 5100 -b 100
```
The results from this script were analzed in Excel to determine the optimal L-value. 


Next, we applied the optimal L-value to define the respective regions using Deliminate_Genome.py:
``` bash
python Deliminate_Genome.py -i Pc2109_SC_OG.csv -l 1100

python Deliminate_Genome.py -i Pc2113_SC_OG.csv -l 1100
```
