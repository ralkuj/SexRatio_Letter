Code related to the letter rebutting some of the conclusions in the article
Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876

Which includes a critique of the conclusions from paper 
Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849

The simulations in first two files (01 and 02) are free of external data. The simulation in last file (04) uses data on family structure from Zietsch et al, data cannot be shared openly.

Files and short description:

01.SNP_simulations_Song&Zhang_20250225.R
Implements the simulation from Song&Zhang to estimate the statistical power to detect a non-zero heritability olf offspring sex-ratio. Implements the original wrong p-value threshold by Song&Zhang  of 5e-8 and a correct threshold of 0.05.

02.SNP_simulations_Song&Zhang_sex-specific_20250525.R
Implements a version of Song&Zhang simulation where parental sex is not random. Uses sample sizes from Zietsch et al to yield power calculations that are trustworthy. First all male-male siblings' offspring simulated, then firstborn to all male-male siblings.

03.SNP_simulations_Plotting_20250225.R
Plots the above simulations.

04.SexRatio_simulations_Normal_liability_20250225.R
Performs a simulation using the relative structure from Zietsch et al to find statistical power for detecting non-null heritability under assumed model with offspring sex allocation drawn from a bernoulli ditribution with probability distributed according to a normal distribution centered aroun 0.5 and varying variance. The variances include the variance used in Song&Zhang, but also smaller and larger variances.
