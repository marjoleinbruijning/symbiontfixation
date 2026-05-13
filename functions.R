
######################################################
## Functions to reproduce figures from:
## The bacterial march to symbiosis: on-ramps and off-ramps
## Bruijning et al.
## Code written by: Marjolein Bruijning
## 2026-05-12
## R-version 4.5.3
######################################################

############################################################################
## Calculate equilibrium frequency given mu, v, s=0
############################################################################
eqprop <- function (mu, v) {
    mu/(mu + v)
}

############################################################################
### Fixation probabilities
############################################################################

## Population genetic model constant population size
fixprob <- function (s,N) {
    (1-exp(-2*s))/(1-exp(-2*s*N))
}

## Fixation probability including host demographic stochasticity
## From: Vindenes et al., 2010
fixprobhostdemo <- function (s,N,sigma,vk,p) {
    (1-exp((-s)/sigma)) / (1-exp((-2*s*N)/sigma))
}

## Fixation probability including fluctuating environments
## From: Cvijovic et al., 2015
fixprobfluct <- function (N,s,sbar,vartau,tau,fluct=FALSE) {
    sstar <- fsstar(s=s,vartau=vartau,tau=tau,N=N,fluct=fluct)
    (2*sbar)/(1-exp(-sbar/sstar))
}

fsstar <- function (s,vartau,tau,N,fluct) {
    xmid <- exp(-s*tau/2)
    tauc <- 1/s

    if (!fluct) {
        (tauc/tau)/(4*N*xmid)
    } else if (fluct) {
      (((s*vartau)^2)/(4*tau)) * (1/(log(N*((s*vartau)^2)*xmid/tauc)))
    }
}


############################################################################
### Calculate demographic variance
## (based on: https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/oik.04708)
############################################################################
demovar <- function (matU,matF) {

    nstages <- ncol(matU)

    ## Stable stage distribution
    w <- matrix(stable.stage(matU+matF),ncol=1) 

    ## reproductive value distribution
    V <- matrix(reproductive.value(matU+matF),ncol=1)
    
    V <- V / sum(t(w) %*% V)
    
    Cu <- Cf <- Ci <- array(NA,dim=c(nstages,nstages,nstages))
    for (i in 1:nstages) {
        p <- matrix(matU[,i],ncol=1)
        mat <- matrix(0,ncol=length(p),nrow=length(p))
        diag(mat) <- p
        Cu[,,i] <- mat - p%*%(t(p))

        Cf[,,i] <- 0
        diag(Cf[,,i]) <- matF[,i]

        Ci[,,i] <- Cu[,,i] + Cf[,,i]
    }

    qvector <- matrix(apply(Ci,3,function(x) t(V) %*% x %*% V),ncol=1)

    t(qvector) %*% w 
}



######################################################
## Simulations for haploid population, with or without flucatuting selection
######################################################

simulateWinner <- function(N, s, fluctuating, times,
                           tau=NA, vartau=NA, sbar=NA) {
    
    if (fluctuating) {

        w1 <- exp(sbar+s)
        w2 <- exp(sbar-s)

        ## Start in one of the two epochs
        if (rbinom(1,1,.5) == 0) {
            signn <- 1
            w <- w1
        } else {
            signn <- -1
            w <- w2
        }

        if (vartau > 0) {
            ## duration of first epoch
            duration <- sample(1:1e4,prob=dnorm(1:1e4,tau,vartau),size=1)
        } else {
            duration <- tau
        }

        ## start at random time in epoch
        nextepoch <- sample.int(round(duration),1)+1

    } else {
        w <- 1+s
    }

    freq <- 1/N

    j <- 1
    continue <- TRUE

    while (continue) {

        if (fluctuating) {
            if (j == nextepoch) {

                ## Switch fitness
                signn <- signn * -1
                if (signn == 1) {
                    w <- w1
                } else if (signn == -1) {
                    w <- w2
                }

                ## When does environment change again?
                if (vartau > 0) {
                    duration <- sample(1:1e4,prob=dnorm(1:1e4,tau,vartau),size=1)
                } else {
                    duration <- round(tau)
                }

                nextepoch <- j+duration
            }

        }

        ## Selection
        psel <- (freq*w) / (freq*w + (1-freq))

        ## Sample next generation
        freq <- rbinom(1,N,psel) / N

        j <- j + 1

        if (freq == 1 | freq == 0 | j == times) {
            continue <- FALSE
        }

    }

    winner <- ifelse(freq == 0,1,2)
    if (j == times) winner <- NA
    if (j == times) warning('Max time reached!')

    return(winner)
}


######################################################
## Simulate dynamics of two competing species with given life history
######################################################
simulateWinnerDemo <- function (N,ss1,ss2,p,matsF,matsP,times) {
    
    ## Initial individuals per stage given stable stage distribution and initial frequencies
    freq <- matrix(c((1-p)*ss1*N,p*ss2*N),nrow=2)

    continue <- TRUE
    i <- 1

    while (continue) {

        ## Projections to next time step
        for (sp in 1:2) { ## for each type
            repr <-  sum(rpois(round(freq[2,sp]),matsF[1,2,sp]))
            freq[,sp] <- matsP[,,sp] %*% freq[,sp,drop=FALSE]
            freq[1,sp] <- freq[1,sp] + repr
        }

        ## Keep total pop size constant
        freq[] <- rmultinom(1,N,freq)

        i <- i + 1
        if (i >  (times-1) | any(colSums(freq) == N)) continue <- FALSE

    }

    winner <- ifelse(colSums(freq)[2] == N,2,1)

    if (i == times) warning('Max time reached!')
    return(winner)
}
