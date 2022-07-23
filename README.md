# CV19BH
Demographic and Viral-Genetic Analyses of COVID-19 Severity in Bahrain Identify Local Risk Factors and a Protective Effect of Polymerase Mutations
1. seq_data_hq_with_mutations_final.csv: all patient data used in analyses including columns for viral SNPs.
2. bahrain_seqs_v5.fasta.gz: genome sequences, names should match Lab.No column in seq_data_hq_with_mutations_final.csv
3. nextclade.tsv: nextclade output after running on bahrain_seqs_v5.fasta.gz
4. problematic_sites_sarsCov2.vcf: obtained from https://github.com/W-L/ProblematicSites_SARS-CoV2.
5. mutation_annotations_082721.csv: counts and consequences of all observed SNPs.
6. Mutation Annotations 061022.ipynb: notebook for generating mutation_annotations_082721.csv.
7. COVID-Mergers.Rmd: notebook used for merging patient information from processed original forms and removing duplicates. Orignal data not provided due to the presence of identifying and private medical information.
8. COVID Vaccination Compare 061022.Rmd: notebook used to perform strain by vaccine analysis.
9. genotype_functions.R: utility functions for various analyses.
10. load_tables.R: process genotypes and filter patient data for presence of COVID-19 severity before analyses.
11. PCA_plots.R: create PCA and other summary plots.
12. nongenetic_anlayses.R: analyses of the effects of nongenetic host factors.
13. individual_mutations.R: association analyses of indvidual viral mutations.
14. rv_burden_plots.R: analysis of burden scores.
15. SFS_plots.R: selection analyses of SNP frequencies. 
