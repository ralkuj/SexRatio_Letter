################################################################################
###
### File: 03.SNP_simulations_Plotting_20250225.R
### Author: Ralf Kuja-Halkola
### Purpose: Plot simulations
### Date: 2024-11-05
### Updated: 2025-02-25
###
################################################################################


############################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849
############################################


############################################
### Repeat of Song&Zhang plot
simulSong_e10m8 <- readRDS(file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_wrong_pval_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSong_e10m8
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/Original_incorrect_pval_20250225.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(0.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Song&Zhang, incorrect p-value threshold',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1  )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### Repeat of Song&Zhang plot, but correct p-value threshold
simulSong <- readRDS(file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_correct_pval_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSong
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/Original_correct_pval_20250225.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Song&Zhang, corrected p-value threshold',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1 )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### More reasonable single-sex model - All offspring
simulSongMalesAll <- readRDS(file='~/Sex ratio/SimulationSNP_male_all_1000repeats_Song_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSongMalesAll
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/MalesOnly_all_20250225.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Males only, all offspring',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1  )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### More reasonable single-sex model - firstborn only
simulSongMalesFirst <- readRDS(file='~/Sex ratio/SimulationSNP_male_firstborn_1000repeats_Song_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSongMalesFirst
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/MalesOnly_firstborn_20250225.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Males only, firstborn',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1 )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
  )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################