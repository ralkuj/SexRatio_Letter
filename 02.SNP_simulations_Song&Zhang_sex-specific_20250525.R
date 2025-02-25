################################################################################
###
### File: 02.SNP_simulations_Song&Zhang_sex-specific_20250525.R
### Author: Ralf Kuja-Halkola
### Purpose: Implement the power simulation by Song & Zhang, but updated with 
###          known sex of individuals whose offspring sex ratio is estimated for
### Date: 2024-11-05
### Updated: 2025-02-25
###
################################################################################


############################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. When theory cannot explain data, the theory needs rethinking. Invited replies to: Orzack SH, Hardy ICW. 2021, and Lehtonen J. 2021. Proc Biol Sci. Mar 31 2021;288(1947):20210304. doi:10.1098/rspb.2021.0304
############################################

############################################
### Libraries used
library(parallel)

### Nr cores used
nrcores <- 50

############################################

################################################################################
### Simulations

############################################
### Based on the simulation from Song&Zhang, updated with known parental sex
### The code is adopted from accompanying code to Song & Zhang, accessed 2025-02-25:
### https://github.com/song88180/Human_sex_ratio/blob/main/02.2_Zietsch_power_analysis.ipynb
### Some comments below kept from Song&Zhang

### Function for calculating cousin pair probabilities for sex combinations
getprobfun <- function(beta=.1,MAF=.1,SR_WT=.52){
  # Calculate sex ratios based on the effect size
  SR_1M <- SR_WT + beta
  SR_2M <- SR_WT + 2*beta
  a <- MAF
  A <- 1-MAF
  # Genotype probabilities
  P_AA <- A^2
  P_Aa <- 2*A*a
  P_aa <- a^2
  # Proportion of the grandparents having genotype XX and XX. This is talking about the shared grandparents of the cousin pair.
  P_GP_AA_AA <- A^4
  P_GP_AA_Aa <- 4*(A^3)*a
  P_GP_AA_aa <- 2*(A^2)*(a^2)
  P_GP_Aa_Aa <- 4*(A^2)*(a^2)
  P_GP_Aa_aa <- 4*A*(a^3)
  P_GP_aa_aa <- a^4
  
  # Proportion of the parents having genotype XX and XX. This is talking about the parents of the cousin pair that are children of the grandparents.
  P_Sib_AA_AA <- P_GP_AA_AA + (0.5^2)*P_GP_AA_Aa + (0.5^4)*P_GP_Aa_Aa
  P_Sib_AA_Aa <- 0.5*P_GP_AA_Aa + (0.5^2)*P_GP_Aa_Aa
  P_Sib_AA_aa <- (0.5^3)*P_GP_Aa_Aa
  P_Sib_Aa_Aa <- (0.5^2)*P_GP_AA_Aa + P_GP_AA_aa + (0.5^2)*P_GP_Aa_Aa + (0.5^2)*P_GP_Aa_aa
  P_Sib_Aa_aa <- (0.5^2)*P_GP_Aa_Aa + 0.5*P_GP_Aa_aa
  P_Sib_aa_aa <- (0.5^4)*P_GP_Aa_Aa + (0.5^2)*P_GP_Aa_aa + P_GP_aa_aa
  
  # Proportion of cousin pairs with different combinations of sexes:
  # Male-Male, Male-Female, Female-Male, Female-Female
  P_Cos_M_M = P_Sib_AA_AA * (SR_WT^2) + P_Sib_AA_Aa * (SR_WT*SR_1M) + 
              P_Sib_AA_aa * (SR_WT*SR_2M) + P_Sib_Aa_Aa * (SR_1M*SR_1M) + 
              P_Sib_Aa_aa * (SR_1M*SR_2M) + P_Sib_aa_aa * (SR_2M*SR_2M)
  
  P_Cos_M_F = P_Sib_AA_AA * (SR_WT*(1-SR_WT)) + 
              0.5 * P_Sib_AA_Aa * ( SR_WT*(1-SR_1M) + (1-SR_WT)*SR_1M ) + 
              0.5 * P_Sib_AA_aa * ( SR_WT*(1-SR_2M) + (1-SR_WT)*SR_2M ) + 
              P_Sib_Aa_Aa * (SR_1M*(1-SR_1M)) + 
              0.5 * P_Sib_Aa_aa * ( SR_1M*(1-SR_2M) + (1-SR_1M)*SR_2M ) + 
              P_Sib_aa_aa * (SR_2M*(1-SR_2M))
  
  P_Cos_F_M = P_Cos_M_F
  
  P_Cos_F_F = P_Sib_AA_AA * (1-SR_WT)^2 + P_Sib_AA_Aa * (1-SR_WT)*(1-SR_1M) + 
              P_Sib_AA_aa * (1-SR_WT)*(1-SR_2M) + P_Sib_Aa_Aa * (1-SR_1M)*(1-SR_1M) + 
              P_Sib_Aa_aa * (1-SR_1M)*(1-SR_2M) + P_Sib_aa_aa * (1-SR_2M)*(1-SR_2M)
  return(c(P_Cos_M_M,P_Cos_M_F,P_Cos_F_M,P_Cos_F_F))
}


### The simulation function
sexratioFun <- function(N=100,beta=.1,MAF=.1,SR_WT=.52-2*.1*.1){
  fortest <- apply( rmultinom( n=N , size=rep(1,4) , prob=getprobfun(beta=beta,MAF=MAF,SR_WT=SR_WT) ) , 1 , sum )
  return(chisq.test(matrix(fortest,2,2))$p.value)
}
sexratioFun(534693,beta=0.097,MAF=0.0025,SR_WT=.514-2*.097*.0025)
sexratioFun(2840847,beta=0.04,MAF=0.5,SR_WT=.514-2*.04*.5)

############################################


############################################
### Single runs

# First born to male siblings, MAF and beta as identified in Song&Zhang. 
# The overall sex ratio of offspring arranged to be as found in Zietsch et al
N <- 534693
n <- 1000
tt0 <- Sys.time()
pvalsSong <- unlist(mclapply( X=rep(N,n) , FUN=sexratioFun , mc.cores=nrcores , 
  beta=.097,MAF=.0025,SR_WT=.514-2*.097*.0025))
tt1 <- Sys.time()
tt1-tt0
mean( pvalsSong < .05 )

# All offspring to male siblings, MAF=.1 and beta=.1
# The overall sex ratio of offspring not handled, and therefore off from true sex ratio
N <- 2840847
n <- 1000
tt0 <- Sys.time()
pvalsSong <- unlist(mclapply( X=rep(N,n) , FUN=sexratioFun , mc.cores=nrcores ,
  beta=0.1,MAF=0.1,SR_WT=.52))
tt1 <- Sys.time()
tt1-tt0
mean( pvalsSong < 5e-8 ) # Incorrect p-value threshold (power=.747 one run)
mean( pvalsSong < 5e-2 ) # Correct p-value threshold  (power=1.000 one run)

############################################


############################################
### Simulation over range of MAF and betas for plotting

### All offspring, overall sex ratio corrected
# Number of offspring pairs to male-male siblings
N <- 2840847
# Number of repeats for each simulation
n <- 1000
MAFseq <- seq(.0025,.1,len=21)
betseq <- seq(.005,0.2,len=21)
# Total number of repeats
n*length(MAFseq)*length(betseq)

# Run simulation
simulSongAll <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simulSongAll[i,j] <- mean( unlist( mclapply( X=rep(N,n) , FUN=sexratioFun , MAF=MAFseq[i] , bet=betseq[j] , SR_WT=.514-2*betseq[j]*MAFseq[i] , mc.cores=nrcores ) ) < 0.05 )
}}
tt1 <- Sys.time()
tt1-tt0
rownames(simulSongAll) <- paste( 'MAF =', MAFseq )
colnames(simulSongAll) <- paste( 'bet =', betseq )

# Save results
saveRDS(object=simulSongAll , file='~/Sex ratio/SimulationSNP_male_all_1000repeats_Song_20241105.Rds' )


### First-born only, overall sex ratio corrected
# Number of parental male-male sibling pairs
N <- 534693
# Number of repeats for each simulation
n <- 1000
MAFseq <- seq(.0025,.1,len=21)
betseq <- seq(.005,0.2,len=21)
# Total number of repeats
n*length(MAFseq)*length(betseq)

# Run simulation
simulSong <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simulSong[i,j] <- mean( unlist( mclapply( X=rep(N,n) , FUN=sexratioFun , MAF=MAFseq[i] , bet=betseq[j] , SR_WT=.514-2*betseq[j]*MAFseq[i] , mc.cores=nrcores ) ) < 0.05 )
}}
tt1 <- Sys.time()
tt1-tt0
rownames(simulSong) <- paste( 'MAF =', MAFseq )
colnames(simulSong) <- paste( 'bet =', betseq )

# Save results
saveRDS(object=simulSong , file='~/Sex ratio/SimulationSNP_male_firstborn_1000repeats_Song_20241105.Rds' )
    

############################################

################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################
