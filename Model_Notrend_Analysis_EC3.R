



#distance model
cat("
    model{
    
    #################
    #Detection model#
    #################
   # priors for distance model; 

pi <- 3.141593

#random year effect for distance sampling model; 
## Priors for hyper-parameters
mu.dd ~ dunif(0,200)
tau.dd <- pow(sigma.dd, -2)
sigma.dd ~ dunif(0, 20)

for(j in 1:n.Sites2){
  eps.dd[j] ~dnorm(0, tau.dd) 
                  }




########################################################
for(t in 1:n.Years){
for(j in 1:n.Sites2){
  sigma[j,t] <- mu.dd + eps.dd[j]
  sigma2[j,t] <- sigma[j,t]*sigma[j,t]
  
  # effective strip width
  esw[j,t] <- sqrt(pi * sigma2[j,t] / 2) 
  p[j,t] <- esw[j,t] / W
  f0[j,t] <- 1/esw[j,t] #assuming all detected on the line
                  }}

########################################################   
for (i in 1:N){ 
  # LIKELIHOOD
  # using zeros trick
  y[i] ~ dunif(0,W) 
  L.f0[i] <- exp(-y[i]*y[i] / (2*sigma2[detectionSite2[i], detectionYear[i]])) * 1/esw[detectionSite2[i], detectionYear[i]] #y are the distances
  nlogL.f0[i] <-  -log(L.f0[i])
  zeros.dist[i] ~ dpois(nlogL.f0[i])
                  }

## Hierarchicla node: esw[j,t]; 
    
########################
#Model of total density#
########################
    
    #priors
    int.d ~ dnorm(0,0.001)    
    
    #random line effect
    line.d.sd ~ dunif(0,10)
    line.d.tau <- pow(line.d.sd,-2)
    for(j in 1:n.Lines){
    random.d.line[j] ~ dnorm(0,line.d.tau)
    }
    

    #random site2 effect
    site2.d.sd ~ dunif(0,10)
    site2.d.tau <- pow(site2.d.sd,-2)
    for(j in 1:n.Sites2){
    random.d.site2[j] ~ dnorm(0,site2.d.tau)
    }

    #random site0 effect
    site0.d.sd ~ dunif(0,10)
    site0.d.tau <- pow(site0.d.sd,-2)
    for(j in 1:n.Sites0){
    random.d.site0[j] ~ dnorm(0,site0.d.tau)
    }
    
    #random time effect
    year.d.sd ~ dunif(0,10)
    year.d.tau <- pow(year.d.sd,-2)
    for(t in 1:n.Years){
    random.d.year[t] ~ dnorm(0,year.d.tau)
    }
    

    #random site2 and time effect
    #s2year.d.sd ~ dunif(0,10)
    #s2year.d.tau <- pow(s2year.d.sd,-2)
    #for(j in 1:n.Sites2){
    #for(t in 1:n.Years){
    #random.d.s2year[j,t] ~ dnorm(0,s2year.d.tau)
    #}
    #}

    #random site0 and time effect
    s0year.d.sd ~ dunif(0,10)
    s0year.d.tau <- pow(s0year.d.sd,-2)
    for(j in 1:n.Sites0){
    for(t in 1:n.Years){
    random.d.s0year[j,t] ~ dnorm(0,s0year.d.tau)
    }
    }
    
    #Model
    for(j in 1:n.Lines){
    for(t in 1:n.Years){

    #linear predictor on density
    log(Density[j,t]) <- int.d + 
    random.d.line[j] + 
    random.d.year[t] + 
    random.d.site0[site0[j]] + 
    random.d.s0year[site0[j],t] +
    random.d.site2[site2[j]]  
    #random.d.s2year[site2[j],t]
    }
    } 
  
    ###################################################
    ###################################################
    ## Linking to number of observed birds; 
        for(j in 1:n.Lines){
        for(t in 1:n.Years){
    
    
        expN_tot_LY[j,t] <- (Density[j,t] * (TransectLength[j,t]/1000 * (W/1000) * 2))
        NuIndivs[j,t] ~ dpois(expNuIndivs[j,t])
        expNuIndivs[j,t] <- Density[j,t] * ((TransectLength[j,t]/1000) * (W/1000) *2) * p[detectionSite2[j],t]

                          }
                          }  
  

    ####################################################
    # N hat - Total individuals in covered area; 

    for (t in 1:n.Years){
    N_tot[t] <- sum(expN_tot_LY[,t])
    N_tot_SD[t] <- N_tot[t]/N_tot[1]
    }


    ####################################################
    # Mean density in (total) covered area 
    
    for (t in 1:n.Years){
    Mean_D[t] <- mean(Density[,t])
    Mean_SD[t] <- Mean_D[t]/Mean_D[1]
    }
    
    ####################################################
    # Mean density in (total) covered area pr. region: 
    
    for(j in 1:n.Sites0){
    for (t in 1:n.Years){
    Mean_DR[j,t] <- sum(Density[,t] * Reg_ind[,j])/sum(Reg_ind[,j])
    #Mean_SDR[j,t] <- Mean_DR[j,t]/Mean_DR[j,1]
    }}

    
    
    ####################################################
    # Mean density in (total) covered area pr. region: 
    
    for(j in 1:n.Sites0){
    EC_index1[j] <- mean(Mean_DR[j,9:13])
    EC_index2[j] <- mean(Mean_DR[j,10:13])
    EC_index3[j] <- mean(Mean_DR[j,11:13])

                        }


    ############################################  
    # Nation wide - based on random effects

    for (t in 1:n.Years){
    Dn[t] <- exp(int.d + random.d.year[t])
    SD[t] <- Dn[t]/Dn[1]  
    }
    
    # Regions
    for(j in 1:n.Sites0){
    for(t in 1:n.Years){
    Dr[j,t] <- exp(int.d + random.d.year[t] + random.d.site0[j] + random.d.s0year[j,t])
    SDr[j, t] <- Dr[j,t]/Dr[j,1]
    }
    }

    ##########################################
    #
    for(j in 1:n.Sites0){
    EC_index_A[j] <- mean(Dr[j,9:13])
    EC_index_B[j] <- mean(Dr[j,10:13])
    EC_index_C[j] <- mean(Dr[j,11:13])

                        }


    }
    ",fill=TRUE,file="NoTrendModel_TrendAnalysis_EC3.txt")
