

######################################################
## Code to reproduce figures from:
## The bacterial march to symbiosis: on-ramps and off-ramps
## Bruijning et al.
## Code written by: Marjolein Bruijning
## 2026-05-12
## R-version 4.5.3
######################################################

## Load all functions
source('functions.R')

## Load required packages
require('fields')

###################################################
## Balance between VT and ET shaping symbiont prevalence
###################################################
n <- 100
mat <- outer(X=seq(0,1,length.out=n),Y=seq(0,1,length.out=n),
             function(X,Y) eqprop(X,Y))

cols <- colorRampPalette(c('#8e0152','#c51b7d','#de77ae',
                           '#f1b6da','#fde0ef','#f7f7f7','#e6f5d0',
                           '#b8e186','#7fbc41','#4d9221','#276419'))(100)

## Plot
fields::image.plot(mat,col=cols,breaks=seq(0,1,l=101),
                   xlab='Symbiont attachment rate',
                   ylab='Symbiont loss rate',
                   cex.axis=1.5,
                   cex.lab=1.5)
lines(0:1,0:1,lwd=2,lty=2)


###################################################
## Equilibrium prevalence for non-commensals
###################################################
par(mfrow=c(1,2),mar=c(2,2,3,1),oma=c(4,1,1,1))

allm <- c(.001,.2,.5) ## transmission rates to show

## Deleterious mutants and strict VT
cols <- c('#fd8d3c','#f03b20','#bd0026')
plot(NA,NA,ylim=c(0,1),xlim=c(-1,0),
      xlab='',
      bty='n',lwd=5,cex.axis=1.2,cex.lab=1.2,
      ylab='',
      yaxt='n',xaxs='i',yaxs='i')
axis(4,labels=NA)
for (i in 1:3) {
    curve(allm[i]/abs(x),add=TRUE,lwd=5,col=cols[i])
}
abline(v=-allm,col=cols,lty=2,lwd=3)
legend('topleft',lwd=3,col=cols,legend=allm,
       title='Environmental\ntransmission',bty='n',cex=.9)

mtext('Strict vertical transmission',3,cex=1.2)


## Beneficial symbionts and without ET
cols <- c('#78c679','#31a354','#006837')
plot(NA,NA,ylim=c(0,1),xlab='',xlim=c(0,1),
     lwd=5,cex.axis=1.2,cex.lab=1.2,ylab='',
     yaxt='n',bty='n',xaxs='i',yaxs='i')
axis(2,cex.axis=1.2)

for (i in 1:3) {
    curve(1-allm[i]/abs(x), add=TRUE,col=cols[i],lwd=5)
}
abline(v=allm,col=cols,lty=2,lwd=3)

legend('bottomright',lwd=3,col=cols,legend=allm,
       title='Loss through\nimperfect vertical\ntransmission',bty='n',cex=.9)
mtext('No environmental transmission',3,cex=1.2)

mtext('Symbiont effect on host fitness',1,
      outer=TRUE,cex=1.2,line=1)
