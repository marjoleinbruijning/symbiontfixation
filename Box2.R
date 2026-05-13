

######################################################
## Code to reproduce figures presented in Box 2, from
## The bacterial march to symbiosis: on-ramps and off-ramps
## Bruijning et al.
## Code written by: Marjolein Bruijning
## 2026-05-12
## R-version 4.5.3
######################################################

## Load all functions
source('functions.R')

## Load required packages
require('parallel')
require('scico')
require('popbio')
require('mcprogress')

## General settings
ncores <- 10
simul <- TRUE ## run all simulations (if FALSE, load from existing files)
plott <- TRUE ## create plots?


######################################################
## Standard pop gen model
######################################################

cat('\r\r \t ------ Classic pop gen ------- \n')

reps <- 5e5
times <- 1e7 # maximum number of time steps


## Settings
allN <- c(100,1e6)
alls <- seq(-.01,.01,l=8)

## run simulations
if (simul) {
    
    df <- expand.grid(N=allN,s=alls)

    usimul <- pmclapply(1:nrow(df),function (k) {

        N <- df$N[k]
        s <- df$s[k]

        fix <- rep(NA,reps)
        for (i in 1:reps) {
            fix[i] <- simulateWinner(N=N,s=s,fluctuating=FALSE,times=times)
        }

        ## How many times has the symbiont reached fixation?
        table(factor(fix,levels=1:2))[2] / length(fix[!is.na(fix)])
        
     }, mc.cores=ncores)

    df$u <- c(unlist(usimul))

    save(df,file='figBox2-a.RData')
    
} else {
    load('figBox2-a.RData')
}


## Create figure
if (plott) {
    cols <- scico(3,palette = 'vanimo')
    x <- seq(-.01,.01,l=100)

    plot(NA,NA,xlim=c(-.01,.01),bty='l',xlab='Endosymbiont effect',
         ylab='Endosymbiont fixation probability',
         cex.lab=1.5,cex.axis=1.5,ylim=c(0,.03))
    abline(v=0,lty=2,lwd=2)

    for (i in 1:length(allN)) {

        ## Add theoretical results
        lines(x,fixprob(x,allN[i]),type='l',lwd=8,col=cols[i])

        ## Add simulation results
        tmp <- df[df$N == allN[i],]
        points(tmp$s,tmp$u,bg=cols[i],lwd=2,cex=2,pch=21,col='grey')
    }

    legend('topleft',col=cols,lwd=5,legend=paste0('N=',allN),
           bty='n',cex=1.5)
}


###################################################
## Fluctuating environments
###################################################

cat('\r\r \t ------ Fluctuating environments ------- \n')

reps <- 1e6
times <- 1e7 # max number of timesteps

## Settings
s <- 10^-2
sbar <- 2e-5
allcv <- c(0,.1)
N <- 1e6

if (simul) {

    ## Create dataframe with all combinations
    n <- 8
    df <- expand.grid(N=N,
                      s=s,
                      tau=round(sqrt(seq(50,3/s,length.out=n))/s),
                      sbar=c(-sbar,sbar),
                      cv=allcv)
    df$vartau <- df$tau * df$cv

    usimul <- pmclapply(1:nrow(df),function (k) {

        vartau <- df$vartau[k]
        tau <- df$tau[k]
        sbar <- df$sbar[k]
        s <- df$s[k]
        N <- df$N[k]

        fix <- rep(NA,reps)
        for (i in 1:reps) {
            fix[i] <- simulateWinner(N=N,s=s,fluctuating=TRUE,times=times,
                                     tau=tau,vartau=vartau,sbar=sbar)
        }

        summ <- table(factor(fix,levels=1:2))
        summ[2] / sum(summ)
        
    }, mc.cores=ncores)

    df$usimul <- c(unlist(usimul))

    save(df,file='figBox2-c.RData')

} else {
    load('figBox2-c.RData')
}

## Plot
if (plott) {
    df2 <- expand.grid(N=N,
                       s=s,
                       tau=seq(min(df$tau),max(df$tau),length.out=100),
                       sbar=c(-sbar,sbar),
                       cv=unique(df$cv))
    df2$vartau <- df2$tau * df2$cv

    cols <- matrix(scico(4, palette = 'hawai'),ncol=2)
    plot(NA,NA,
         xlim=(df2$s[1]*c(min(df2$tau),max(df2$tau)))^2,
         ylim=c(10,1000),
         log='y',
         xlab=bquote("Scaled timescale of fluctuations [" * (tau %*% s)^2 * "]"),
         ylab='Scaled endosymbiont fixation probability (N x u)',
         bty='l',
         cex.axis=1.5,cex.lab=1.5)
    abline(h=2*sbar*N,lwd=2,lty=2)


    for (i in 1:length(unique(df$sbar))) {
        for (j in 1:length(unique(df$cv))) {
            
            inc2 <- df2$sbar == unique(df$sbar)[i] & df2$cv == unique(df$cv)[j]
            pred <- fixprobfluct(
                sbar=df2$sbar[inc2],
                s=df2$s[inc2],
                tau=df2$tau[inc2],
                N=df2$N[inc2],
                vartau=df2$vartau[inc2],
                fluct=ifelse(unique(df$cv)[j] == 0, FALSE, TRUE))
            
            lines((df2$tau[inc2]*df2$s[inc2])^2,pred*df2$N[inc2],col=cols[j,i],lwd=8)

            ## Add simulations
            inc <- df$sbar == unique(df$sbar)[i] & df$cv == unique(df$cv)[j]
            points((df$tau[inc]*df$s[inc])^2,(df$usimul[inc])*df$N[inc],bg=cols[j,i],pch=21,
                   cex=2,col='grey',lwd=2)
        }
    }

    legend('topleft',col=cols,lwd=3,
           legend=c(bquote(bar(s) == .(-sbar) ~ ", " ~ CV[tau] == .(allcv[1])),
                    bquote(bar(s) == .(-sbar) ~ ", " ~ CV[tau] == .(allcv[2])),
                    bquote(bar(s) == .(sbar) ~ ", " ~ CV[tau] == .(allcv[1])),
                    bquote(bar(s) == .(sbar) ~ ", " ~ CV[tau] == .(allcv[2]))),
           bty='n',cex=1.5)
}


####################################################
## Host demographic stochasticity
####################################################

cat('\r\r \t ------ Host life history ------- \n')

reps <- 5e5
times <- 1e7 # maximum number of time steps

## Settings
N <- 1e6 ## pop size
alls <- c(.0001,.001,.01) ## symbiont effects
p <- 1/N ## initial freq symbiont carrying hosts

if (simul) {

    ## Create dataframe with combinations of growth and survival
    n <- 3
    scenarios <- expand.grid(surv=(seq((0.1),(0.5),length.out=n)),
                             growth=(seq((.1),(0.9),length.out=n)),
                             s=alls)
    
    ## Set fecundity so that lambda=1
    scenarios$f <- -((scenarios$growth-1)*scenarios$surv^2+
                     (2-scenarios$growth)*scenarios$surv-1)/(scenarios$growth*scenarios$surv) 

    scenarios$usimul <- scenarios$demvar <- NA

    usimul <- pmclapply(1:nrow(scenarios), function(sc) {

        fitnesseffect <- 1+scenarios$s[sc]

        g1 <- scenarios$growth[sc]
        s1 <- scenarios$surv[sc]
        f1 <- scenarios$f[sc]

        ## Symbiont free hosts
        ## Fill matrices
        matF1 <- matP1 <- matF2 <- matP2 <- matrix(0,nrow=2,ncol=2)
        diag(matP1) <- c(s1*(1-g1),s1)
        matP1[2,1] <- g1*s1
        matF1[1,2] <- f1

        ## Symbiont carrying host
        g2 <- g1
        s2 <- s1 * fitnesseffect
        f2 <- f1 * fitnesseffect
        diag(matP2) <- c(s2*(1-g2),s2)
        matP2[2,1] <- g2*s2
        matF2[1,2] <- f2

        ## Combine both species into one array
        matsF <- matsP <- array(NA,dim=c(2,2,2))
        matsP[,,1] <- matP1
        matsP[,,2] <- matP2
        matsF[,,1] <- matF1
        matsF[,,2] <- matF2

        ## Start with stable distributions
        ss1 <- stable.stage(matP1+matF1)
        ss2 <- stable.stage(matP2+matF2)
        ss1 <- ss1 / sum(ss1)
        ss2 <- ss2 / sum(ss2)
        
        demvar <- c(demovar(matF=matsF[,,1],matU=matsP[,,1]))

        fix <- rep(NA,reps)
        for (r in 1:reps) {
            fix[r] <- simulateWinnerDemo(N=N,ss1=ss1,ss2=ss2,
                                         p=p,matsF=matsF,
                                         matsP=matsP,times=times)
        }

        summ <- table(factor(fix,levels=1:2))
        prop <- summ[2] / sum(summ)
        c(demvar,prop)
        
    }, mc.cores=ncores)

    scenarios$usimul <- sapply(usimul,function(x) x[2])
    scenarios$demvar <- sapply(usimul,function(x) x[1])

    save(scenarios,file='figBox2-b.RData')
    
} else {
    load('figBox2-b.RData')
}


if (plott) {

    cols <- scico(length(alls),palette='glasgow')
    plot(NA,NA,xlim=c(.5,10),
         lwd=5,bty='l',xlab='Host demographic variance',
         ylab='Scaled endosymbiont fixation probability (N x u)',
         cex.lab=1.5,cex.axis=1.5,ylim=c(1,1e5),log='xy')

    for (i in 1:length(alls)) {
        curve(fixprobhostdemo(s=alls[i],N=N,
                              sigma=x,vk=1,p=1/N)*N,add=TRUE,
              lwd=8,col=cols[i])
        abline(h=2*alls[i]*N,col=cols[i],lwd=2,lty=2)

        inc <- scenarios$s == alls[i]
        points(scenarios$demvar[inc],scenarios$usimul[inc]*N,bg=cols[i],
               pch=21,col='grey',lwd=2,cex=2)
    }

    legend('bottomleft',bty='n',col=cols,
           legend=paste0('s=',alls),lwd=5,cex=1.5)
}


########################################################################################
