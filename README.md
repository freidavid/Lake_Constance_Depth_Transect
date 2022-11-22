## Introgression from extinct species facilitates adaptation to its vacated niche


All scripts to generate the results (starting with processing of raw sequence data, mapping, calculation of genotype likelihoods and the analyses done) can be found in the folder lsf_submit_scripts. Scripts were run according to the sequence described in the file on the ETHZ cluster Euler.

All scripts to visualize and finalize the data are in the folder R_scripts.


### lsf_submit_scripts.


#### fastp.lsf
Used to trim poly-G tails from raw reads.


#### Seqprep.lsf
Used to merge overlapping forward and reverse reads.


#### BWA.lsf
Used to map the data to the Alpine whitefish reference genome (De-Kayne et al. 2020).


#### Processing_bam.lsf
Used to process the mapped bam files.


#### Angsd_per_chromosom.lsf
Used to calculate genotypelikelihoods (per chromosome).


#### merging_single_beagle.lsf
Used to merge the single beagle files of each chromosomes into one genome-wide file with genotype likelihoods.


#### PCAngsd.lsf
Used to do a genomic PCA, admixture analysis and selection scan.


#### d_stats.sh
This script submits all the submit script d-stat._final.lsf to do all needed sinlge d-stat calculations.


#### maf_at_positions.lsf
Used to calculated allele frequencies at given positions (here the positions that are significant in the selection scan).


#### maf_genomewide.lsf
Used to calcualate all genome wide allele frequencies.


#### get_overlapping_genes.lsf
Used to extract the genes that are overlapping with the positions that are signficant in the selection scan.







