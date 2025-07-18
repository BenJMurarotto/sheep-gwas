🧬 GWAS Analysis of Unknown Sheep Phenotype
This project explores the genetic basis of an unknown phenotype in sheep using Genome-Wide Association Studies (GWAS). Conducted as part of the GENE552 course at the University of New England, the analysis used genotypic data (48,570 SNPs) from 300 sheep and corresponding residual phenotype values (adjusted for environmental effects).

🔬 Methods
Data Preprocessing: Genotype matrix centered by allele frequency; individuals kept as rows.

Statistical Model: Linear models fitted SNP-by-SNP using lm() in R.

Visualisations: Manhattan and Q-Q plots were generated to identify significant associations.

📈 Key Results
One SNP on chromosome 16 (OAR16_57327865.1) exceeded the Bonferroni threshold (p = 7.33e-07).

Moderate genomic inflation observed (GIF = 1.408).

🧠 Candidate Genes
Several genes near the significant SNP may influence traits such as:

ARHGEF39: Cell growth and motility.

TESK1: Cytoskeletal dynamics and reproductive development.

CA9: pH regulation and hypoxia response.

CD72: B-cell signaling and immune regulation.
