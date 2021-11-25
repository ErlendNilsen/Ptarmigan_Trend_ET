

#### EXAMPLE - RUNNING THE HDS MODELS BASED ON GBIF-DATA; 

load("Bugs_GBIF_data2009_2019.RData")


###########################################################
#### First example; running model M_random on the total population 
#### GBIF-data 2009-2019; 


# MCMC settings
niter <- 30000
nthin <- 1
nburn <- 20000
nchains <- 3    

# Setting parameters to monitor            

#params <- c("b.df.0")

params <- c("mu.dd", "int.d","line.d.sd","year.d.sd",
            "site.d.sd","site2.d.sd", "random.d.year",
            "Mean_D", "Mean_SD", "DN", "SD", "Dr", "SDr",  
            "N_tot", "N_tot_SD", "beta", "Lam_SD")

inits1 <- function() {list(int.d=rnorm(1, mean=log(5), sd=0.1), 
                           mu.dd=runif(1, 80, 90))}

out_GBIF_NoTrend_Total <- jagsUI::jags(bugs.data_GBIF_Total, inits=inits1, params, 
                        model.file="NoTrendModel_TrendAnalysis2020.txt",
                        n.chain=nchains, n.iter=niter, 
                        n.burnin=nburn, parallel = TRUE, 
                        DIC=FALSE, codaOnly=c("Deviance", "Density"))


###########################################################
#### Second example; running model M_dynamic on the adult population 
#### GBIF-data 2009-2019; 


# MCMC settings
niter <- 300
nthin <- 1
nburn <- 200
nchains <- 3    

# Setting parameters to monitor            

#params <- c("b.df.0")

params <- c("r_mean", "N_tot", "N_tot_SD", "Mean_D", "Mean_SD", "r_s0", "r_T")

inits1 <- function() {list(r_mean=runif(1, 0.03, 0.07),
                           mu.dd=runif(1, 80, 90))}

out_GBIF_Trend_Adults <- jagsUI::jags(bugs.data_GBIF_Adults, inits=inits1, params, 
                                       model.file="Model_Trend_Analysis2020.txt",
                                       n.chain=nchains, n.iter=niter, 
                                       n.burnin=nburn, parallel = TRUE, 
                                       DIC=FALSE, codaOnly=c("Deviance", "Density"))


