library(msm)
library(RColorBrewer)
library(emdbook)

#Function to calculate 95% quantiles
cred <- function(x){
  return(quantile(x, c(0.025, 0.975), na.rm=TRUE))
  
#Function to calculate 50% quantiles
}
cred.50 <- function(x){
  return(quantile(x, c(0.25, 0.75)), na.rm=TRUE)
}
#calculate 95% quantiles for columns
cred.2 <- function(x){
  return(apply( x, 2, cred), na.rm=TRUE)
}
#calculate 50% quantiles for columns
cred.3 <- function(x){
  return(apply( x, 2, cred.50), na.rm=TRUE)
}
#Function to calculate the medians of a column
colMedian <- function(x){
  return(apply(x, 2, median), na.rm=TRUE)
}

#recruitment function --> Beverton-Holt or Ricker, parameterized for adult to smolt
rec_func <- function(S=7500, alpha=10.8, beta=0.032, sigma_R=0.0, s_juv=1, form='B-H'){
  #calculate recruitment deviation using mean juvenile survival and lognormal bias correction
  epsilon <- exp((rnorm(1, log(s_juv)-((sigma_R^2)/2), sigma_R)))
  if(form=='B-H'){
    #calculate recruits (Beverton-Holt)
    R <- ((alpha*S)/(1+beta*S))*epsilon
  }
  if(form=='Ricker'){
    #calculate recruits (Ricker)
    R <- alpha*S*exp(-beta*S)
  }
  return(R)
}
rec_func()

#mating function --> model frequency-dependent reproductive success of jacks and hooknose males
mate_func <- function(jacks=10, hooknose=100, form='linear', slope=0.5, eq=0.1, steep=0.99, plot=TRUE){
  #predicted ratio of jack to hooknose reproductive success
  pred_mate_ratio <- c()
  ratio_vec <- seq(from=0, to=0.5, by=0.01)
  if(form=='linear'){
    intercept <- eq*(1-slope)
    if(jacks < 1){
      pred_mate_ratio <- 0
    }
    if (jacks >= 1){
      pred_mate_ratio <- intercept + slope*(jacks/hooknose)
    }
    if(plot==TRUE){
      plot(0, type='n', ylim=c(0, max(ratio_vec)), xlim=c(0,max(ratio_vec)), xlab=NA, ylab=NA)
      for(i in 1:length(ratio_vec)){
        points(ratio_vec[i], intercept + slope*(ratio_vec[i]), pch=16)
        abline(0,1)
        mtext(side=1, 'Jack:Hooknose ratio in spawning population', line=2.75, cex=2)
        mtext(side=2, 'Jack:Hooknose reproductive success', line=2.5, cex=2)
      }
    }
  }
  if(form=='saturating'){
    alpha <-  ((1-steep)/(4*steep*eq))*eq
    beta <-  (5*steep-1)/(4*steep*eq)
    pred_mate_ratio <- (jacks/hooknose)/(alpha + beta*(jacks/hooknose))
    if(plot==TRUE){
      plot(0, type='n', ylim=c(0, max(ratio_vec)), xlim=c(0,max(ratio_vec)))
      for(i in 1:length(ratio_vec)){
        points(ratio_vec[i], (ratio_vec[i])/(alpha + beta*(ratio_vec[i])))
        abline(0,1)
      }
    }    
  }
  return(pred_mate_ratio)
}
mate_func()

#Growth function --> von Bertallanfy growth function based on previous size
sVBGM_s<-function(size,Linf_vB=750,k_vB=0.3,sd_vB=0.2) { 
  return(size+k_vB*(Linf_vB-size)*exp(rnorm(length(Linf_vB),-0.5*sd_vB^2,sd_vB)))
} 

## probabilistic maturation reaction norm based on size and age
PMRN <-function(age=5,size=100,a=500,b=-50,sd=50) {
  return(1/(1+exp(-(size-(a+b*age))/sd)))
}

#Harvest selectivity function
HARV<-function(size,a=a_harv,sd=sd_harv){ 
  exp(-(log(size)-log(a))^2/(2*sd^2)) 
}

## adult survival function (D)
adult_surv <- function(input, s_adult=0.8){
  #calculate survival as a bernouli trial for each individual with a probability equal to the survival term
  x <- rbinom(n=nrow(input), size=1, prob=s_adult)
  #save the row indices for individuals that survived
  out <- matrix(input[which(x==1),], byrow=F, ncol=4)
  return(out)
}

#Plot all functions
plot_fun <- function(S_R_alpha=60, S_R_beta=0.00017, Linf_vB= 700, k_vB=0.3, a_PMRN=800, b_PMRN=-100, sd_PMRN=25, eq=0.1, Linf.vs.k=-0.0005, Linf.vs.PMRN=0.3, sel_alpha=475, sel_sd=1){
  pdf(file='func_plot.pdf')
  par(mfrow=c(3,3), oma=c(10,10,10,10), mar=c(2,2,2,2))
  #Recruitment
  s_vec <- c(1:100000)
  r_vec <- vector(length=length(s_vec))
  for(i in 1:length(s_vec)){
    r_vec[i] <- rec_func(S=s_vec[i], alpha=S_R_alpha, beta=S_R_beta, sigma_R=0.0, s_juv=1, form='B-H')
  }
  plot(s_vec, r_vec, type='l', xlab='Female spawners', ylab='Smolts', lwd=2)
  mtext(side=1,'Female spawners', cex=0.75, line=2)
  mtext(side=2,'Smolts', cex=0.75, line=2)
  legend("bottomright", "A", bty="n") 
  #growth
  growth_vec <- vector(length=4)
  growth_vec[1] <- 100
  for(i in 2:length(growth_vec)){
    growth_vec[i] <- sVBGM_s(growth_vec[i-1] ,Linf_vB=Linf_vB,k_vB=k_vB,sd_vB=0.0)
  }
  plot(1:4, growth_vec[c(1:4)], type='o', pch=16, col='black', lwd=2, xlab='Age', ylab='Length (mm)', xaxt='n')
  axis(side=1, at=c(1,2,3,4), labels=c('0', '1', '2', '3'))
  mtext(side=1,'Ocean age', cex=0.75, line=2)
  mtext(side=2,'Length (mm)', cex=0.75, line=2)
  legend("bottomright", "B", bty="n") 
  
  #Maturation probabilities
  length_vec <- c(100:900)
  p_mat <- matrix(nrow=length(length_vec), ncol=3)
  for(j in 1:ncol(p_mat)){
    for(i in 1:length(length_vec)){
      p_mat[i,j] <- PMRN(age=j+2,size=length_vec[i],a=a_PMRN,b=b_PMRN,sd=sd_PMRN)
    }
  }
  plot(length_vec, p_mat[,1], type='l', col='gray80', ylim=c(0,1), xlim=c(100, 700), lwd=2, xlab='Length (mm)', ylab='Maturation probability')
  abline(v=length_vec[which(p_mat[,1]==0.5)], col='gray80', lty=2)
  points(length_vec[which(p_mat[,1]==0.5)], 0.5, col='gray80')
  points(length_vec, p_mat[,2], type='l', col='gray50', lwd=2)
  abline(v=length_vec[which(p_mat[,2]==0.5)], col='gray50', lty=2)
  points(length_vec[which(p_mat[,2]==0.5)], 0.5, col='gray50')
  points(length_vec, p_mat[,3], type='l', col='black', lwd=2)
  abline(v=length_vec[which(p_mat[,3]==0.5)], col='black', lty=2)
  points(length_vec[which(p_mat[,3]==0.5)], 0.5, col='black')
  mtext(side=1,'Length (mm)', cex=0.75, line=2)
  mtext(side=2,'Maturation probability', cex=0.75, line=2)
  legend("bottomright", "C", bty="n") 
  
  #Growth parameter tradeoff
  L_inf_vec <- c(400,1000)
  k_vec <- vector(length=length(L_inf_vec))
  for(i in 1:length(k_vec)){
    k_vec[i] <- k_vB +  Linf.vs.k*(L_inf_vec[i]-Linf_vB)
  }
  
  plot(L_inf_vec, k_vec, type='l', lwd=2, ylab='Growth rate coefficient', xlab='Asymptotic size', ylim=c(0.1,0.5), xlim=c(400,1000), xaxs='i')
  abline(v=Linf_vB, lwd=0.1)
  abline(h=k_vB, lwd=0.1)
  mtext(side=1,'Asymptotic length (mm)', cex=0.75, line=2)
  mtext(side=2,'Growth rate coefficient', cex=0.75, line=2)
  legend("topright", "D", bty="n") 
  
  #PMRN asymptotic length tradeoff
  PMRN_vec <- c(400:1200)
  L_inf_vec <- vector(length=length(PMRN_vec))
  
  for(i in 1:length(L_inf_vec)){
    L_inf_vec[i] <- Linf_vB + Linf.vs.PMRN*(PMRN_vec[i]-a_PMRN)
  }
  plot(PMRN_vec, L_inf_vec, type='l', lwd=2, ylab='Asymptotic size', xlab='PMRN intercept', xlim=c(400,1200), ylim=c(400, 1000), xaxs='i')
  abline(h=Linf_vB, lwd=0.1)
  abline(v=a_PMRN, lwd=0.1)
  mtext(side=1,'PMRN intercept', cex=0.75, line=2)
  mtext(side=2,'Asymptotic size', cex=0.75, line=2)
  legend("bottomright", "E", bty="n") 
  
  #Harvest selectivity
  #Gillnet selectivity function
  size_vec <- seq(100,1000,1)
  sel_sd_vec <- lseq(from=0.1, to=2, length.out=10)
  col.vec <- gray.colors(n=length(sel_sd_vec))
  sel_vec <- matrix(nrow=length(size_vec), ncol=length(sel_sd_vec))
  plot(NA, xlim=c(250,550), type='l', lwd=2, xlab='size(mm)', ylab='Harvest selectivity', ylim=c(0,1))
  for(i in 1:length(size_vec)){
    for(j in length(sel_sd_vec):1){
      sel_vec[i,j] <-  HARV(size_vec[i], a=sel_alpha, sd=sel_sd_vec[j])
      points(size_vec, sel_vec[,j], type='l', lwd=2, lty=1, col=col.vec[j])
    }
  }
  mtext(side=1,'Length (mm)', cex=0.75, line=2)
  mtext(side=2,'Harvest selectivity', cex=0.75, line=2)
  legend("bottomright", "F", bty="n") 
  
  #Frequency dependent mating
  jack_vec <- c(1:500)
  mate_vec <- matrix(ncol=10, nrow=length(jack_vec))
  slope_vec <- seq(0, 0.9, 0.1)
  col.vec <- gray.colors(n=10)
  
  for(j in 1:length(slope_vec)){
    for(i in 1:nrow(mate_vec)){
      mate_vec[i,j] <- mate_func(jacks=jack_vec[i], hooknose=1000, form='linear', slope=slope_vec[j], eq=eq, steep=0.99, plot=FALSE)
    }
  }
  plot(NA, xlim=c(0,0.5), ylim=c(0,0.5), xaxs='i', ylab='Jack:hooknose fertilizations', xlab='Jack:hooknose')
  for(j in 1:length(slope_vec)){
    points(jack_vec/1000, mate_vec[,j], type='l', lwd=2, col=col.vec[j])
  }
  abline(0,1, lty=2, lwd=1)
  abline(v=0.1, lty=3, lwd=1)
  mtext(side=1,'Jack:hooknose present', cex=0.75, line=2)
  mtext(side=2,'Jack:hooknose fertilizations', cex=0.75, line=2)
  legend("bottomright", "G", bty="n") 
  
  dev.off()
}
plot_fun()

rec_vec <- c(0.4, 0.6, 0.8)
slope_vec <- c(0.1, 0.3, 0.5)
growth_vec <- c(0.2)
int_vec <- c(0.01)
harvest_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
sel_vec <- rev(lseq(from=0.1, to=2, length.out=10))

#List for storing the total jack prevalence by brood year
jack_mat <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the prevalence of jacks by return year
jack_mat_return <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the prevalence of jacks by brood year among males only
jack_mat_male <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the prevalence of jacks by return year among males only
jack_mat_return_male <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the harvest
catch_mat <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the realized exploitation rate
u_real_mat <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))
#List for storing the harvest, exclusing jacks
non_jack_catch_mat <- replicate(length(rec_vec), list(replicate(length(slope_vec), list(matrix(ncol=length(harvest_vec), nrow=length(sel_vec))))))

set.seed(666)
  for(x in 1:length(sel_vec)){
    for(h in 1:length(harvest_vec)){
      for(w in 1:length(slope_vec)){
        for(r in 1:length(rec_vec)){
          for(g in 1:length(growth_vec)){
            for(k in 1:length(int_vec)){
              burnin <- 1   # Deterministic burning
              n_year <- 1000   # Total length of simulation
              n_age <- 5      # Maximum number of age classes
              n_sim <- 10      # Number of stochastic simulations to run
              n_sex <- 2      # Number of sexes (2)
              n_sire <- 2     # Number of sire types (2)
              Linf<- 700			# Reference value for asymptotic length
              k_vb<- 0.3     # Reference value for growth rate
              Linf.vs.k <- -0.0005		# Slope of correlation of Linf vs k
              Linf.vs.PMRN <- 0.3     # Slope of correlation of Linf vs PMRN intercept
              sd.g_g<- growth_vec[g]	# Growth variation
              SaS<-100		              # Smolt size
              s_adult <- c(1, 0.8, 0.8) #Age specific survival probabilities
              s_juv <- 0.2    # Early marine survival    
              b_mat<- -100	  # PMRN slope
              sd_mat <- 25    # Width of the PMRN
              S_R_alpha <- 60     # Productivity for S-R (Coho, Barrowman et al 2003)
              S_R_beta <-0.00017  # Capacity for S-R (Coho, Barrowman et al 2003)
              sigma_R_g <- rec_vec[r]    # Recruitment variation                
              sex_ratio <- 0.5     # Sex ratio mean
              sex_var_g <- 0.0     # Variance in the sex ratio
              PMRN_male <- 800     # Reference value for the intercept of the PMRN (male)
              sd_male <- 200       # Standard deviation for the intercept of the PMRN (male)
              male_lower <- 200    # Lower bound for the intercept of the PMRN (male)
              male_upper <- 2000   # Upper bound for the intercept of the PMRN (male)
              PMRN_female <- 800   # Reference value for the intercept of the PMRN (female)
              sd_female <- 200     # Standard deviation for the intercept of the PMRN (female)
              female_lower <- 200  # Lower bound for the intercept of the PMRN
              female_upper <- 2000 # Upper bound for the intercept of the PMRN
              mate_eq <-int_vec[k] # Equilibrium jack:hooknose mating frequency
              prop_jack_sire_init <- mate_eq   #Initial proportion of jack-sired individuals
              lin_mate_slope <- slope_vec[w]   #Slope of the linear frequency-dependent mating function
              sat_mate_steep <- NA # Steepness of the saturating linear frequency-dependent mating function
              recomb_error <- 0.071 #Mutation, segregation, recombination kernel
              plots<-TRUE #Diagnostic plots
              male_exp <- 1 #Allometric exponent relating body size to matine success for males
              female_exp <- 1.77 #Allometric exponent relating body size to fecundity for females
              #Harvest parameters
              harvest_rate <- harvest_vec[h] #Fully selected harvest rate
              harvest <-TRUE #Is harvest ocurring or not
              harvest_burn <- 500 #How many additional years after the burnin until harvest is allowed?
              sel_alpha <- 475 #Center of the selectivity curve
              sel_sd <- sel_vec[x] #Width of the selectivity curve
              
              #loop over number of stochastic iterations
              #create lists for mature and immature individuals
              imm_list <- replicate(n_year,list(replicate(n_age,list(replicate(n_sex,list(replicate(n_sire, list(replicate(n_sim,list())))))))))
              mat_list <- replicate(n_year,list(replicate(n_age,list(replicate(n_sex,list(replicate(n_sire, list(replicate(n_sim,list())))))))))
              
              #lists for storing harvested and escaped individuals
              mat_list_catch <- replicate(n_year,list(replicate(n_age,list(replicate(n_sex,list(replicate(n_sire, list(replicate(n_sim,list())))))))))
              mat_list_esc <- replicate(n_year,list(replicate(n_age,list(replicate(n_sex,list(replicate(n_sire, list(replicate(n_sim,list())))))))))
              
              #matrix for storing the realized exploitation rates
              u_real <- matrix(nrow=n_year, ncol=n_sim)
              
              #matrix for storing proportional reproductive success of jacks over time
              prop_jack_sire_mat <- matrix(nrow=n_year, ncol=n_sim)
              
              for(l in 1:n_sim){
                #loop over years
                for (y in 1:n_year){
                  #year 1
                  if(y==1){
                    #initiate rectuits based on the capcity assumed in the S-R function
                    recruits <- S_R_alpha/S_R_beta 
                    #Calculate the sex ratio of the recruits as 0.5 plus some small error
                    prop_female <- rtnorm(1, sex_ratio, 0, lower=0, upper=1)
                    
                    #loop through ages
                    for(a in 1:n_age){
                      #first age class:age 1 individuals (fry/recruits)
                      if(a==1){
                        #determine female recruits, hooknose sire
                        imm_list[[y]][[a]][[1]][[1]][[l]] <- matrix(ncol=4, nrow=round(recruits*prop_female*(1-prop_jack_sire_init)))
                        #determine male recruits, hooknose sire
                        imm_list[[y]][[a]][[2]][[1]][[l]] <- matrix(ncol=4, nrow=round(recruits*(1-prop_female)*(1-prop_jack_sire_init)))
                        #determine female recruits, jack sire 
                        imm_list[[y]][[a]][[1]][[2]][[l]] <- matrix(ncol=4, nrow=round(recruits*prop_female*prop_jack_sire_init))
                        #determine male recruits, jack sire
                        imm_list[[y]][[a]][[2]][[2]][[l]] <- matrix(ncol=4, nrow=round(recruits*(1-prop_female)*prop_jack_sire_init))
                        
                        #each individual is monitored for four traits, which are stored in a four column matrix:
                        # 1) phnotypic intercept of the PMRN (inherited from parents) --> column 1
                        # 2) Linfinity drawn from population wide-distribution --> column 2
                        # 3) Brody growth coeffient (deterministic function of Linfinity) --> column 3
                        # 4) length (calculated via stochastic von-bertalanfy) --> column 4
                        
                        #loop through sex 
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #if individuals are male and were sired by a jack, they inherit a PMRN intercept for jacks
                            if(s==1){
                              imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(length(imm_list[[y]][[a]][[s]][[j]][[l]][,1]), mean=PMRN_female, sd=sd_female, lower=female_lower, upper=female_upper)
                              imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]- PMRN_female)
                              #otherwise they inherit a standard value
                            }
                            if(s==2){
                              imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(length(imm_list[[y]][[a]][[s]][[j]][[l]][,1]), mean=PMRN_male, sd=sd_male, lower=male_lower, upper=male_upper)
                              imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]-PMRN_male)
                            }
                            
                            #calculate Brody growth coefficients as a deterministic function of Linifinity
                            imm_list[[y]][[a]][[s]][[j]][[l]][,3] <- k_vb+Linf.vs.k*(imm_list[[y]][[a]][[s]][[j]][[l]][,2]-Linf)
                            #Calculate length as function of stochastic von-Bertalanfy equation and initial size
                            imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- SaS
                          }
                        }
                      }
                      #age 2 individuals (smolts--> assuming everyone spends only 1 year in the lake after gravel year)
                      if(a > 1 & a < 3){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #First calculate survival --> for year 1 to initiate the simulation, assume that the number age 2 individuals is equal to the number of age 1 individuals in the same year times the survival term
                            #To monitor individuals, survival is implemented by randomly sampling individuals with a probability according to the adult survival term 
                            imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y]][[a-1]][[s]][[j]][[l]], s_adult=1)
                            #Individual growth of surviving individuals is calculated according to the stochastic von-bertalanfy equation (for year 1, previous size is based on the size of the previous age class in the same year)
                            #imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=0)
                            imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- SaS
                            
                            }
                        }
                      }
                      #Age 3 and older individuals (assume individuals of these ages are in the ocean)
                      if(a > 2 & a < n_age){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #only perform calculations if there are any immature individuals of the previous age class remaining in the population
                            if(is.null(dim(imm_list[[y]][[a-1]][[s]][[j]][[l]]))==FALSE){
                              # calculate surviving individuals
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y]][[a-1]][[s]][[j]][[l]], s_adult=s_adult[a-2])
                              #calculate growth of surviving individuals according to stochastic von-Bertalanfy eq.
                              imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=0)
                              #calculate maturation probabilities according to PMRN, and determine maturation as a bernouli distributed random variable with this probability
                              if(a ==3 & s==1){
                                mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=0.01)
                              }else{
                                mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=PMRN(age=a,size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],a=imm_list[[y]][[a]][[s]][[j]][[l]][,1],b=b_mat,sd=sd_mat))
                              }
                              #'mature' is a boolean identifying which individuals in the age class matured and which did not
                              # use these indices to identify the rows containing mature versus immature individuals and split into mature and immature lists accordingly
                              mat_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][which(mature==1),], ncol=4, byrow=F)
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][-which(mature==1),], ncol=4, byrow=F)
                            }
                          }
                        }
                      }
                      #maximum age class (plus group) --> any individuals of this age mature with 100% probability
                      if(a==n_age){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #only perform calculations if there are any immature individuals of the previous age class remaining in the population
                            if(is.null(dim(imm_list[[y]][[a-1]][[s]][[j]][[l]]))==FALSE){
                              # calculate surviving individuals
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y]][[a-1]][[s]][[j]][[l]], s_adult=s_adult[a-2])
                              #calculate growth of surviving individuals according to stochastic von-Bertalanfy eq.
                              imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=0)
                              #calculate maturation probabilities according to PMRN, and determine maturation as a bernouli distributed random variable with this probability
                              mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=1)
                              #'mature' is a boolean identifying which individuals in the age class matured and which did not
                              # use these indices to identify the rows containing mature versus immature individuals and split into mature and immature lists accordingly 
                              mat_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][which(mature==1),], ncol=4, byrow=F)
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][-which(mature==1),], ncol=4, byrow=F)
                            }
                          }    
                        }
                      }
                    }
                  }   
                  if(y > 1){#for years greater than one
                    if(y <= burnin){ #Keep variance terms set to zero within the burnin
                      sigma_R <- 0
                      sd.g <- 0
                      sex_var <- 0
                    }
                    if(y > burnin){ # Beyond the burnin, use the actual variance terms defined
                      sigma_R <- sigma_R_g
                      sd.g <- sd.g_g
                      sex_var <- sex_var_g
                    }
                    
                    total_return <- 0  #Calculate the total return (mature individuals), initiate at 0
                    for (a in 3:n_age){
                      for (s in 1:n_sex){
                        for (j in 1:n_sire){
                          if(is.null(dim(mat_list[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_return <- total_return + nrow(mat_list[[y-1]][[a]][[s]][[j]][[l]])
                          }
                        }
                      }
                    }
                    
                    if(harvest == TRUE & y > (burnin+harvest_burn)){ #additional harvest burnin
                      #Individuals are harvested 
                      n_harvest <- total_return*harvest_rate #number determined by the return size and the harvest rate
                      F_catch <- -log(1-harvest_rate) #calculate fishing mortality rate from the supplied target exploitation rate
                      for (a in 3:n_age){
                        for (s in 1:n_sex){
                          for (j in 1:n_sire){
                            if(is.null(dim(mat_list[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              #probability of harvest determined by target fishing mortality and individual size
                                catch <- rbinom(n=nrow(mat_list[[y-1]][[a]][[s]][[j]][[l]]), size=1, prob=(1-exp(-F_catch*HARV(size=mat_list[[y-1]][[a]][[s]][[j]][[l]][,4], a=sel_alpha, sd=sel_sd))))
                              #'catch' is a boolean identifying which individuals in the age class are caught
                              # use these indices to identify the rows containing harvested versus escaped individuals and split into mature and immature lists accordingly 
                              mat_list_catch[[y-1]][[a]][[s]][[j]][[l]] <- matrix(mat_list[[y-1]][[a]][[s]][[j]][[l]][which(catch==1),], ncol=4, byrow=F)
                              mat_list_esc[[y-1]][[a]][[s]][[j]][[l]] <-  matrix(mat_list[[y-1]][[a]][[s]][[j]][[l]][-which(catch==1),], ncol=4, byrow=F)

                            }
                          }
                        }
                      }
                      
                      #calculate the total harvest, initiate at 0
                      total_harvest <- 0
                      for (a in 3:n_age){
                        for (s in 1:n_sex){
                          for (j in 1:n_sire){
                            if(is.null(dim(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_harvest <- total_harvest + nrow(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]])
                            }
                          }
                        }
                      }
                      #calculate the total escapement, initiate at 0
                      total_esc <- 0
                      for (a in 3:n_age){
                        for (s in 1:n_sex){
                          for (j in 1:n_sire){
                            if(is.null(dim(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_esc <- total_esc + nrow(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]])
                            }
                          }
                        }
                      }
                      
                      #Calculate true exploitation rate
                      u_real[y,l] <- total_harvest/total_return
                    
                      
                     #Calculate the total escapement for males and females
                      total_esc_male <- 0
                      total_esc_female <- 0
                      
                      #Calculate the total harvest for males and females, initiate at 0
                      total_catch_male <- 0
                      total_catch_female <- 0
                      
                      for (a in 3:n_age){
                        for (s in 1:n_sex){
                          for (j in 1:n_sire){
                            if(s==1){
                              if(is.null(dim(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                                total_esc_female <- total_esc_female + nrow(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]) #looping through all ages, sex, and sire types, add the female fish that matured and were not caught in the previous year to the total female escapement
                              }
                              if(is.null(dim(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                                total_catch_female <- total_catch_female + nrow(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]) #looping through all ages, sex, and sire types, add the female fish that matured and were caught in the previous year to the total female harvest
                              }
                            }
                            if(s==2){
                              if(is.null(dim(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                                total_esc_male <- total_esc_male + nrow(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]) #looping through all ages, sex, and sire types, add the male fish that matured and were not caught in the previous year to the total male escapement
                              }
                              if(is.null(dim(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                                total_catch_male <- total_catch_male + nrow(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]) #looping through all ages, sex, and sire types, add the male fish that matured and were caught in the previous year to the total male harvest
                              }
                            }
                          }
                        }
                      }
                      #replace NaNs in the abundance with 0s
                      total_esc_female[is.na(total_esc_female)] <- 0
                      total_catch_female[is.na(total_catch_female)] <- 0
                      
                      total_esc_male[is.na(total_esc_male)] <- 0
                      total_catch_male[is.na(total_catch_male)] <- 0
                      
                      #create empty vectors for the age composition
                      age_comp_esc_male <- vector(length=n_age)
                      age_comp_esc_female <- vector(length=n_age)
                      
                      age_comp_catch_male <- vector(length=n_age)
                      age_comp_catch_female <- vector(length=n_age)
                      
                      #Loop through age and sex to calculate the age composition of the harvest and escapement (for both sire types)
                      for(a in 1:n_age){
                        for(s in 1:n_sex){
                          if(s==2){
                            if(is.null(dim(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              age_comp_esc_male[a] <- (nrow(mat_list_esc[[y-1]][[a]][[s]][[1]][[l]]) + nrow(mat_list_esc[[y-1]][[a]][[s]][[2]][[l]]))/total_esc_male
                            }
                            if(is.null(dim(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              age_comp_catch_male[a] <- (nrow(mat_list_catch[[y-1]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y-1]][[a]][[s]][[2]][[l]]))/total_catch_male
                            }
                          }
                          if(s==1){
                            if(is.null(dim(mat_list_esc[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              age_comp_esc_female[a] <- (nrow(mat_list_esc[[y-1]][[a]][[s]][[1]][[l]]) + nrow(mat_list_esc[[y-1]][[a]][[s]][[2]][[l]]))/total_esc_female
                            }
                            if(is.null(dim(mat_list_catch[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                              age_comp_catch_female[a] <- (nrow(mat_list_catch[[y-1]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y-1]][[a]][[s]][[2]][[l]]))/total_catch_female
                            }
                          }
                        }
                      }
                      
                      #replace NaNs in the age composition with 0s
                      age_comp_esc_female[is.na(age_comp_esc_female)] <- 0
                      age_comp_esc_male[is.na(age_comp_esc_male)] <- 0
                      
                      age_comp_catch_female[is.na(age_comp_catch_female)] <- 0
                      age_comp_catch_male[is.na(age_comp_catch_male)] <- 0
                      
                      #calculate average annual age-specific harvest selectivity on males (use for calculating fitness ratio later)
                      sel_a <- vector(length=(n_age-2))
                      for (a in 3:n_age){
                        if(is.null(dim(mat_list[[y-1]][[a]][[s]][[j]][[l]]))==FALSE){
                          if(nrow(mat_list[[y-1]][[a]][[s]][[j]][[l]]) > 0){
                            #For each age class of males, calculate the number that were harvested divided by the total that reached maturity
                            sel_a[a-2] <- (nrow(mat_list_catch[[y-1]][[a]][[2]][[1]][[l]]) + nrow(mat_list_catch[[y-1]][[a]][[2]][[2]][[l]]))/(nrow(mat_list[[y-1]][[a]][[2]][[1]][[l]]) + nrow(mat_list[[y-1]][[a]][[2]][[2]][[l]]))
                          }
                        }
                        else{
                          sel_a[a-2] <- 0
                        }
                      }
                    }
                    #If there is no harvest or too early in the simulation...
                    if(harvest == FALSE | y <= burnin+harvest_burn){
                      mat_list_esc <- mat_list
                      sel_a <- rep(0, length=(n_age-2))
                    }
                    #create object to store the spawning population, partitioned by jacks, hooknoses and females
                    #jacks --> merge jacks that were sired by jacks and those sired by hooknoses
                    jacks <- rbind(mat_list_esc[[y-1]][[3]][[2]][[1]][[l]], mat_list_esc[[y-1]][[3]][[2]][[2]][[l]])
                    
                    #hooknose males --> initiate at age 4, merge across sire types
                    #include a vector for calculating the age composition of hooknose males, used later fitness calculations
                    hook_age_comp <- vector(length=(n_age-3))
                    hook <- matrix(rbind(mat_list_esc[[y-1]][[4]][[2]][[1]][[l]], mat_list_esc[[y-1]][[4]][[2]][[2]][[l]]), ncol=4)
                    hook_age_comp[1] <- nrow(matrix(rbind(mat_list_esc[[y-1]][[4]][[2]][[1]][[l]], mat_list_esc[[y-1]][[4]][[2]][[2]][[l]]), ncol=4))
                
                    
                    if(n_age > 4){
                      #loop over remaining ages
                      for (a in 5:n_age){
                        hook <- matrix(rbind(hook, mat_list_esc[[y-1]][[a]][[2]][[1]][[l]], mat_list_esc[[y-1]][[a]][[2]][[2]][[l]]), ncol=4) #Grab age 5 males of both sire types
                        hook_age_comp[a-3] <- nrow(mat_list_esc[[y-1]][[a]][[2]][[1]][[l]]) + nrow(mat_list_esc[[y-1]][[a]][[2]][[2]][[l]]) #Grab the number of age 5 males of both sire types
                      }
                    }

                    #check
                    sum(hook_age_comp)==nrow(hook) #(TRUE?)
                    #calculate hooknose age composition 
                    hook_age_comp <- hook_age_comp/nrow(hook)
                    
                    #Create an object for storing females --> initiate at age 3 (should be rare but allowed), merge across sire types
                    females <- matrix(rbind(mat_list_esc[[y-1]][[3]][[1]][[1]][[l]], mat_list_esc[[y-1]][[3]][[1]][[2]][[l]]), ncol=4)
                    females <- cbind(females, rep(3, nrow(females)))
                    for (a in 4:n_age){
                      #Grab all females from the mature list, add a column identifying the age of each female
                      females <- matrix(rbind(females, cbind(mat_list_esc[[y-1]][[a]][[1]][[1]][[l]], rep(a, nrow(mat_list_esc[[y-1]][[a]][[1]][[1]][[l]]))), cbind(mat_list_esc[[y-1]][[a]][[1]][[2]][[l]], rep(a, nrow(mat_list_esc[[y-1]][[a]][[1]][[2]][[l]])))), ncol=5)
                    }
                    
                    #Calculate recruits in a given year as a function of mature females and S-R parameters
                    recruits <- rec_func(S=nrow(females), alpha=S_R_alpha, beta=S_R_beta, sigma_R=sigma_R, s_juv=s_juv, form='B-H')
                    
                    #calculate the proportional reproductive output of each female as a function of length
                    fecundity <- rbinom(nrow(females), round(recruits),  ((females[,4]^female_exp)/sum(females[,4]^female_exp)))

                    #Calculate the proportion of recruits sired by jacks according to negatively frequency-dependent mating function
                    pred_mate <- mate_func(jacks=nrow(jacks), hooknose=nrow(hook), form='linear', eq=mate_eq, slope=lin_mate_slope, steep=sat_mate_steep, plot=FALSE)
                    if(dim(hook)[1] > 0){
                      prop_jack_sire <- pred_mate/(pred_mate+1) #Assuming there is at least 1 hooknose male, calculate proportion from the ratio returned by the FDS functions
                    }
                    if(dim(hook)[1] == 0){
                      prop_jack_sire <- 1.0 # If there are no hooknose males, 100% of offspring are sired by jacks
                    }
                    
                    #Save the propotion sired by jacks to a matrix
                    prop_jack_sire_mat[y,l] <- prop_jack_sire
                    
                    #calculate individual hooknoses probability of mating according to their body size
                    hook_prob <- ((hook[,4]^male_exp)/sum(hook[,4]^male_exp))*(1- prop_jack_sire_mat[y,l])
                    
                    #calculate individual jacks probability of spawning
                    jack_prob <- rep(prop_jack_sire_mat[y,l]/nrow(jacks), nrow(jacks))
                    
                    #combine jacks and hooknoses + their mating probabilities into single objects for all males
                    males <- rbind(jacks, hook)
                    
                    #Log transform male mating probabilities
                    male_prob <- c(log(jack_prob), log(hook_prob))
                    
                    #calculate sex ratio
                    prop_female <- rtnorm(1, sex_ratio, sex_var, lower=0, upper=1)
                    
                    #assign recruits to sexes and sires according to each male's probability of mating --> males identifued and selected by row number
                    sire <- sample((1:nrow(males)), size=nrow(females), prob=exp(male_prob), replace=TRUE)
                    
                    # create matrix with each female (column 1), the row number of the male who sired her offspring (column 2), the number of offspring produced (column 3)
                    #and a blank column (4) to identify whether or not the father was a jack
                    parents <- matrix(c((1:nrow(females)), sire, fecundity, rep(NA, nrow(females))), nrow=nrow(females), ncol=4, byrow=F)
                    
                    #Specify column 4 as a binary code identifying whether the father is a jack or not
                    parents[parents[,2]<= nrow(jacks),4] <- 2 #2 = jack father
                    parents[parents[,2]>= nrow(jacks),4] <- 1 #1 = hooknose father
                    
                    #expand the parents matrix so each row represents an individual offspring with identifiers for the row number of the mother and father, and identifying whether the father was a jack or not
                    if(nrow(parents>0)){
                      parents_2 <- parents[rep(1:nrow(parents), parents[,3]),-3]
                      parents_2 <-cbind(parents_2, rbinom(nrow(parents_2), 1, prop_female)+1)
                      
                    }else{
                      parents_2 <- parents
                    }
                    
                    #calculate weighted probability of survival for hooknose males (includes harvest)
                    hook_surv <- vector(length=length(hook_age_comp))
                    for(a in 1:length(hook_age_comp)){
                      hook_surv[a] <- hook_age_comp[a]*(s_adult[a+1]^a)*(1-sel_a[a+1])
                    }
                    hook_surv <- sum(hook_surv)
                    
                    #calculate fitness ratio
                    jack_mate <-prop_jack_sire/nrow(jacks)
                    hook_mate <- (1-prop_jack_sire)/nrow(hook)

                    #loop through ages
                    for(a in 1:n_age){
                      #first age class:age 1 individuals (fry/recruits)
                      if(a==1){
                        #determine female recruits, hooknose sire
                        imm_list[[y]][[a]][[1]][[1]][[l]] <- matrix(ncol=4, nrow=nrow(as.matrix(parents_2[parents_2[,3]==1 & parents_2[,4]==1,])))
                        
                        #determine male recruits, hooknose sire
                        imm_list[[y]][[a]][[2]][[1]][[l]] <- matrix(ncol=4, nrow=nrow(as.matrix(parents_2[parents_2[,3]==1 & parents_2[,4]==2,])))
                        
                        #determine female recruits, jack sire 
                        imm_list[[y]][[a]][[1]][[2]][[l]] <- matrix(ncol=4, nrow=nrow(as.matrix(parents_2[parents_2[,3]==2 & parents_2[,4]==1,])))
                        
                        #determine male recruits, jack sire
                        imm_list[[y]][[a]][[2]][[2]][[l]] <- matrix(ncol=4, nrow=nrow(as.matrix(parents_2[parents_2[,3]==2 & parents_2[,4]==2,])))
                        
                        #each individual is monitored for four traits, which are stored in a four column matrix:
                        # 1) phnotypic intercept of the PMRN (inherited from parents) --> column 1
                        # 2) Linfinity drawn from population wide-distribution --> column 2
                        # 3) Brody growth coeffient (deterministic function of Linfinity) --> column 3
                        # 4) length (calculated via stochastic von-bertalanfy) --> column 4
                        #loop through sex and sire identity
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #if individuals are male and were sired by a jack, they inherit a value from their father
                            if(j==2 & s==2){
                              if(length(jacks)>0){
                                imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), males[matrix(parents_2[parents_2[,3]==2 & parents_2[,4]==2,], ncol=4)[,2],1], recomb_error*PMRN_male , lower=male_lower, upper=male_upper)
                                imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]- PMRN_male)
                                
                              }
                            }
                            #males sired by hooknose males inherit values from their father
                            if(j==1 & s==2){
                              if(length(hook)>0){
                                imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), males[matrix(parents_2[parents_2[,3]==1 & parents_2[,4]==2,],ncol=4)[,2],1], recomb_error*PMRN_male, lower=male_lower, upper=male_upper)
                                imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]- PMRN_male)
                                
                              } 
                            }
                            #females of both sire types inherit a value from their mothers
                            if(j==1 & s==1){
                              if(length(females)>0){
                                imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), females[matrix(parents_2[parents_2[,3]==1 & parents_2[,4]==1,],ncol=4)[,1],1], recomb_error*PMRN_female, lower=female_lower, upper=female_upper)
                                imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]- PMRN_female)
                              }
                            }
                            if(j==2 & s==1){
                              if(length(females)>0){
                                imm_list[[y]][[a]][[s]][[j]][[l]][,1] <- rtnorm(nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), females[matrix(parents_2[parents_2[,3]==2 & parents_2[,4]==1,],ncol=4)[,1],1], recomb_error*PMRN_female, lower=female_lower, upper=female_upper)
                                imm_list[[y]][[a]][[s]][[j]][[l]][,2] <- Linf + Linf.vs.PMRN*(imm_list[[y]][[a]][[s]][[j]][[l]][,1]- PMRN_female)
                              }
                            }
                            
                            #calculate growth coefficients as a deterministic function of Linifinity
                            imm_list[[y]][[a]][[s]][[j]][[l]][,3] <- k_vb+Linf.vs.k*(imm_list[[y]][[a]][[s]][[j]][[l]][,2]-Linf)
                            #Initialize length
                            imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- SaS
                          }
                        }
                      }
                      #age 2 individuals
                      if(a > 1 & a < 3){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #First calculate survival --> assume survival is 100%, S-R function integrates mortality through the first ocean age for computational efficiency
                            imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y-1]][[a-1]][[s]][[j]][[l]], s_adult=1)
                            #Individual growth of surviving individuals is calculated according to the stochastic von-bertalanfy equation
                            imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=sd.g)
                            imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- SaS
                          }
                        }
                      }
                      #Age 3 and older individuals (assume individuals of these ages are in the ocean)
                      if(a > 2 & a < n_age){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #only perform calculations if there are any immature individuals of the previous age class remaining in the population
                            if(is.null(dim(imm_list[[y-1]][[a-1]][[s]][[j]][[l]]))==FALSE){
                              # calculate surviving individuals
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y-1]][[a-1]][[s]][[j]][[l]], s_adult=s_adult[a-2])
                              #calculate growth of surviving individuals according to stochastic von-Bertalanfy eq.
                              imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=sd.g)
                              #calculate maturation probabilities according to PMRN, and determine maturation as a bernouli distributed random variable with this probability
                              if(a == 3 & s==1){
                                mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=0.01)
                              }else{
                                mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=PMRN(age=a,size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],a=imm_list[[y]][[a]][[s]][[j]][[l]][,1],b=b_mat,sd=sd_mat))
                              }                      
                              #'mature' is a boolean identifying which individuals in the age class matured and which did not
                              # use these indices to identify the rows containing mature versus immature individuals and split into mature and immature lists accordingly
                              mat_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][which(mature==1),], ncol=4, byrow=F)
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][-which(mature==1),], ncol=4, byrow=F)
                            }
                          }
                        }
                      }
                      #maximum age class (plus group) --> any individuals of this age mature with 100% probability
                      if(a==n_age){
                        #loop through sex and sire
                        for(s in 1:n_sex){
                          for(j in 1:n_sire){
                            #only perform calculations if there are any immature individuals of the previous age class remaining in the population
                            if(is.null(dim(imm_list[[y-1]][[a-1]][[s]][[j]][[l]]))==FALSE){
                              # calculate surviving individuals
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- adult_surv(input=imm_list[[y-1]][[a-1]][[s]][[j]][[l]], s_adult=s_adult[a-2])
                              #calculate growth of surviving individuals according to stochastic von-Bertalanfy eq.
                              imm_list[[y]][[a]][[s]][[j]][[l]][,4] <- sVBGM_s(size=imm_list[[y]][[a]][[s]][[j]][[l]][,4],Linf_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,2],k_vB=imm_list[[y]][[a]][[s]][[j]][[l]][,3],sd_vB=sd.g)
                              #calculate maturation probabilities according to PMRN, and determine maturation as a bernouli distributed random variable with this probability
                              mature <- rbinom(n=nrow(imm_list[[y]][[a]][[s]][[j]][[l]]), size=1, prob=1)
                              #'mature' is a boolean identifying which individuals in the age class matured and which did not
                              # use these indices to identify the rows containing mature versus immature individuals and split into mature and immature lists accordingly 
                              mat_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][which(mature==1),], ncol=4, byrow=F)
                              imm_list[[y]][[a]][[s]][[j]][[l]] <- matrix(imm_list[[y]][[a]][[s]][[j]][[l]][-which(mature==1),], ncol=4, byrow=F)
                            }
                          }    
                        }
                      }
                    }
                  }
                }
              }

              if(plots==TRUE){
                pdf(file=paste('sigma_R=',sigma_R_g,'sigma_g=',sd.g_g ,'slope=',slope_vec[w] ,'harvest=',harvest,'n_age=',n_age,'mate_eq=',mate_eq,'h_rate=',harvest_rate,'sel=',sel_sd,'.pdf'))
                par(mfrow=c(4,4), mar=c(2.5,2.5,2.5,2.5), oma=c(1,1,1,1))
                
                #plot output
                #jack abundance in the total population (total and separated by sire type)
                #jacks sired by hooknose males
                jack_spawn_1 <- matrix(nrow=n_year, ncol=n_sim)
                #jacks sired by jacks
                jack_spawn_2 <- matrix(nrow=n_year, ncol=n_sim)
                #total jacks
                jack_spawn_total <- matrix(nrow=n_year, ncol=n_sim)
                
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    jack_spawn_1[y,l] <- nrow(mat_list[[y]][[3]][[2]][[1]][[l]])
                    jack_spawn_2[y,l] <- nrow(mat_list[[y]][[3]][[2]][[2]][[l]])
                    jack_spawn_total[y,l] <- nrow(mat_list[[y]][[3]][[2]][[1]][[l]]) + nrow(mat_list[[y]][[3]][[2]][[2]][[l]])
                  }
                }
                #Plot time series of jack abundance across sire types and total, mean and 95% quantiles across simulations
                #total jack abundance
                plot(apply(jack_spawn_total, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(jack_spawn_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Jack abundance', main='Jack abundance (return year)', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Abundance', cex=0.85, line=2)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(jack_spawn_total, 1, cred)[1,], rev(apply(jack_spawn_total, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by jacks
                points(apply(jack_spawn_2, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(jack_spawn_2, 1, cred)[1,], rev(apply(jack_spawn_2, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by non-jacks
                points(apply(jack_spawn_1, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(jack_spawn_1, 1, cred)[1,], rev(apply(jack_spawn_1, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,100, max=250, alpha=90), border=FALSE)
                legend(x=100, y=max(apply(jack_spawn_total, 1, cred)[2,]), legend=c('jack-sired', 'hooknose-sired', 'total jacks'), cex=0.65, col=c('red', 'blue', 'black'), lty=1, bty='n')
                
                
                #jack proportion in the spawning population (total and separated by sire type)
                #jacks sired by hooknose males
                jack_p_spawn_1 <- matrix(nrow=n_year, ncol=n_sim)
                #jacks sired by jacks
                jack_p_spawn_2 <- matrix(nrow=n_year, ncol=n_sim)
                #total jacks
                jack_p_spawn_total <- matrix(nrow=n_year, ncol=n_sim)
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                total_spawners <- matrix(nrow=n_year, ncol=n_sim)
                for(y in 1:(n_year-1)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_spawners[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_spawners[y,l] <- total_spawners[y,l] + nrow(mat_list_esc[[y]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                    jack_p_spawn_1[y,l] <- nrow(mat_list_esc[[y]][[3]][[2]][[1]][[l]])/total_spawners[y,l]
                    jack_p_spawn_2[y,l] <- nrow(mat_list_esc[[y]][[3]][[2]][[2]][[l]])/total_spawners[y,l]
                    jack_p_spawn_total[y,l] <- (nrow(mat_list_esc[[y]][[3]][[2]][[1]][[l]]) + nrow(mat_list_esc[[y]][[3]][[2]][[2]][[l]]))/total_spawners[y,l]
                  }
                }
                
                #replace NaNs in the jack proportions with 0s
                jack_p_spawn_1[is.na(jack_p_spawn_1)] <- 0
                jack_p_spawn_2[is.na(jack_p_spawn_2)] <- 0
                jack_p_spawn_total[is.na(jack_p_spawn_total)] <- 0
                
                #Plot time series of jack proportion across sire types and total, mean and 95% quantiles across simulations
                #total jack abundance
                plot(apply(jack_p_spawn_total, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(jack_p_spawn_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Jack proportion', main='Jack prop in escapement', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(jack_p_spawn_total, 1, cred)[1,], rev(apply(jack_p_spawn_total, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by jacks
                points(apply(jack_p_spawn_2, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(jack_p_spawn_2, 1, cred)[1,], rev(apply(jack_p_spawn_2, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by non-jacks
                points(apply(jack_p_spawn_1, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(jack_p_spawn_1, 1, cred)[1,], rev(apply(jack_p_spawn_1, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,100, max=250, alpha=90), border=FALSE)
                legend(x=100, y=max(apply(jack_p_spawn_total, 1, cred)[2,]), legend=c('jack-sired', 'hooknose-sired', 'total jacks'), cex=0.65, col=c('red', 'blue', 'black'), lty=1, bty='n')
                
                #jack proportion in the harvest (total and separated by sire type)
                #jacks sired by hooknose males
                jack_p_catch_1 <- matrix(nrow=n_year, ncol=n_sim)
                #jacks sired by jacks
                jack_p_catch_2 <- matrix(nrow=n_year, ncol=n_sim)
                #total jacks
                jack_p_catch_total <- matrix(nrow=n_year, ncol=n_sim)
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total catching population in a given year
                total_catch <- matrix(nrow=n_year, ncol=n_sim)
                non_jack_catch <- matrix(nrow=n_year, ncol=n_sim)
                
                for(y in (burnin+1):(n_year-1)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_catch[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_catch[y,l] <- total_catch[y,l] + nrow(mat_list_catch[[y]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                    if(is.null(dim(mat_list_catch[[y]][[3]][[2]][[1]][[l]]))==FALSE){
                      jack_p_catch_1[y,l] <- nrow(mat_list_catch[[y]][[3]][[2]][[1]][[l]])/total_catch[y,l]
                    }else{
                      jack_p_catch_1[y,l] <- 0
                    }
                    if(is.null(dim(mat_list_catch[[y]][[3]][[2]][[2]][[l]]))==FALSE){
                      jack_p_catch_2[y,l] <- nrow(mat_list_catch[[y]][[3]][[2]][[2]][[l]])/total_catch[y,l]
                    }else{
                      jack_p_catch_2[y,l] <- 0
                    }
                    if(is.null(dim(mat_list_catch[[y]][[3]][[2]][[1]][[l]]))==FALSE & is.null(dim(mat_list_catch[[y]][[3]][[2]][[2]][[l]]))==FALSE){
                      jack_p_catch_total[y,l] <- (nrow(mat_list_catch[[y]][[3]][[2]][[1]][[l]]) + nrow(mat_list_catch[[y]][[3]][[2]][[2]][[l]]))/total_catch[y,l]
                    }else{
                      jack_p_catch_total[y,l] <- 0
                    }
                  }
                }
                
                non_jack_catch[y,l] <- total_catch[y,l] - jack_p_catch_total[y,l]*total_catch[y,l]
                
                #replace NaNs in the jack proportions with 0s
                jack_p_catch_1[is.na(jack_p_catch_1)] <- 0
                jack_p_catch_2[is.na(jack_p_catch_2)] <- 0
                jack_p_catch_total[is.na(jack_p_catch_total)] <- 0
                
                #Plot time series of jack proportion across sire types and total, mean and 95% quantiles across simulations
                #total jack abundance
                plot(apply(jack_p_catch_total, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(jack_p_catch_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Jack proportion', main='Jack prop in harvest', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(jack_p_catch_total, 1, cred)[1,], rev(apply(jack_p_catch_total, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by jacks
                points(apply(jack_p_catch_2, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(jack_p_catch_2, 1, cred)[1,], rev(apply(jack_p_catch_2, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by non-jacks
                points(apply(jack_p_catch_1, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(jack_p_catch_1, 1, cred)[1,], rev(apply(jack_p_catch_1, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,100, max=250, alpha=90), border=FALSE)
                legend(x=100, y=max(apply(jack_p_catch_total, 1, cred)[2,]), legend=c('jack-sired', 'hooknose-sired', 'total jacks'), cex=0.65, col=c('red', 'blue', 'black'), lty=1, bty='n')
                
                #jack proportion in the total population (total and separated by sire type) by return year
                #jacks sired by hooknose males
                jack_p_run_1 <- matrix(nrow=n_year, ncol=n_sim)
                #jacks sired by jacks
                jack_p_run_2 <- matrix(nrow=n_year, ncol=n_sim)
                #total jacks
                jack_p_run_total <- matrix(nrow=n_year, ncol=n_sim)
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total runing population in a given year
                total_run <- matrix(nrow=n_year, ncol=n_sim)
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_run[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_run[y,l] <- total_run[y,l] + nrow(mat_list[[y]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                    jack_p_run_1[y,l] <- nrow(mat_list[[y]][[3]][[2]][[1]][[l]])/total_run[y,l]
                    jack_p_run_2[y,l] <- nrow(mat_list[[y]][[3]][[2]][[2]][[l]])/total_run[y,l]
                    jack_p_run_total[y,l] <- (nrow(mat_list[[y]][[3]][[2]][[1]][[l]]) + nrow(mat_list[[y]][[3]][[2]][[2]][[l]]))/total_run[y,l]
                  }
                }
                
                #replace NaNs in the jack proportions with 0s
                jack_p_run_1[is.na(jack_p_run_1)] <- 0
                jack_p_run_2[is.na(jack_p_run_2)] <- 0
                jack_p_run_total[is.na(jack_p_run_total)] <- 0
                
                
                #Now just the jack proprtion among males
                #jack proportion in the total population (total and separated by sire type)
                #jacks sired by hooknose males
                jack_p_run_1_male <- matrix(nrow=n_year, ncol=n_sim)
                #jacks sired by jacks
                jack_p_run_2_male <- matrix(nrow=n_year, ncol=n_sim)
                #total jacks
                jack_p_run_total_male <- matrix(nrow=n_year, ncol=n_sim)
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total runing population in a given year
                total_run_male <- matrix(nrow=n_year, ncol=n_sim)
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_run_male[y,l] <- 0
                    for(a in 1:n_age){
                      for(j in 1:n_sire){
                        #conditional to avoid errors from empty lists
                        if(is.null(dim(mat_list[[y]][[a]][[2]][[j]][[l]]))==FALSE){
                          total_run_male[y,l] <- total_run_male[y,l] + nrow(mat_list[[y]][[a]][[2]][[j]][[l]]) 
                        }
                      }
                    }
                    jack_p_run_1_male[y,l] <- nrow(mat_list[[y]][[3]][[2]][[1]][[l]])/total_run_male[y,l]
                    jack_p_run_2_male[y,l] <- nrow(mat_list[[y]][[3]][[2]][[2]][[l]])/total_run_male[y,l]
                    jack_p_run_total_male[y,l] <- (nrow(mat_list[[y]][[3]][[2]][[1]][[l]]) + nrow(mat_list[[y]][[3]][[2]][[2]][[l]]))/total_run_male[y,l]
                  }
                }
                
                #replace NaNs in the jack proportions with 0s
                jack_p_run_1_male[is.na(jack_p_run_1)] <- 0
                jack_p_run_2_male[is.na(jack_p_run_2)] <- 0
                jack_p_run_total_male[is.na(jack_p_run_total)] <- 0
                
                #Plot time series of jack proportion across sire types and total, mean and 95% quantiles across simulations
                #total jack abundance
                plot(apply(jack_p_run_total, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(jack_p_run_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Jack proportion', main= 'Total jack prop (return year)', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(jack_p_run_total, 1, cred)[1,], rev(apply(jack_p_run_total, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by jacks
                points(apply(jack_p_run_2, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(jack_p_run_2, 1, cred)[1,], rev(apply(jack_p_run_2, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by non-jacks
                points(apply(jack_p_run_1, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(jack_p_run_1, 1, cred)[1,], rev(apply(jack_p_run_1, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,100, max=250, alpha=90), border=FALSE)
                legend(x=100, y=max(apply(jack_p_run_total, 1, cred)[2,]), legend=c('jack-sired', 'hooknose-sired', 'total jacks'), cex=0.65, col=c('red', 'blue', 'black'), lty=1, bty='n')
                
                
                #Look at spawners and recruits over time
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                total_spawners <- matrix(nrow=n_year, ncol=n_sim)
                total_recruits <- matrix(nrow=n_year, ncol=n_sim)
                total_catch <- matrix(nrow=n_year, ncol=n_sim)
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_spawners[y,l] <- 0
                    total_recruits[y,l] <- 0
                    total_catch[y,l] <- 0
                    for(s in 1:n_sex){
                      for(j in 1:n_sire){
                        if(is.null(dim(imm_list[[y]][[1]][[s]][[j]][[l]]))==FALSE){
                          total_recruits[y,l] <- total_recruits[y,l] + nrow(imm_list[[y]][[1]][[s]][[j]][[l]]) 
                        }
                      }
                    }
                    
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_spawners[y,l] <- total_spawners[y,l] + nrow(mat_list_esc[[y]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_catch[y,l] <- total_catch[y,l] + nrow(mat_list_catch[[y]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                  }
                }
                
                #Plot time series of spawners and recruits
                plot(apply(total_spawners, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(total_recruits, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Abundance', main='Spawners and recruits', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Abundance', cex=0.85, line=2)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(total_spawners, 1, cred)[1,], rev(apply(total_spawners, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                
                points(apply(total_recruits, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(total_recruits, 1, cred)[1,], rev(apply(total_recruits, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                
                points(apply(total_catch, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(total_catch, 1, cred)[1,], rev(apply(total_catch, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                legend(x=100, y=1.1*max(apply(total_recruits[-1,], 1, cred)[2,]), legend=c('spawners', 'recruits', 'harvest'), cex=0.65, col=c('black', 'red', 'blue'), lty=1, bty='n')
                
                #Jack proportion by brood year
                n_brood <- n_year-n_age
                #jack proportion in the spawning population (total and separated by sire type)
                #jacks sired by hooknose males
                jack_p_brood_1 <- matrix(nrow=n_brood, ncol=n_sim)
                #jacks sired by jacks
                jack_p_brood_2 <- matrix(nrow=n_brood, ncol=n_sim)
                #total jacks
                jack_p_brood_total <- matrix(nrow=n_brood, ncol=n_sim)
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                total_return <- matrix(nrow=n_brood, ncol=n_sim)
                for(y in 1:(n_brood)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_return[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(is.null(dim(mat_list[[y+a]][[a]][[s]][[j]][[l]]))==FALSE){
                            total_return[y,l] <- total_return[y,l] + nrow(mat_list[[y+a]][[a]][[s]][[j]][[l]]) 
                          }
                        }
                      }
                    }
                    jack_p_brood_1[y,l] <- nrow(mat_list[[y+3]][[3]][[2]][[1]][[l]])/total_return[y,l]
                    jack_p_brood_2[y,l] <- nrow(mat_list[[y+3]][[3]][[2]][[2]][[l]])/total_return[y,l]
                    jack_p_brood_total[y,l] <- (nrow(mat_list[[y+3]][[3]][[2]][[1]][[l]]) + nrow(mat_list[[y+3]][[3]][[2]][[2]][[l]]))/total_return[y,l]
                  }
                }
                
                #replace NaNs in the jack proportions with 0s
                jack_p_brood_1[is.na(jack_p_brood_1)] <- 0
                jack_p_brood_2[is.na(jack_p_brood_2)] <- 0
                jack_p_brood_total[is.na(jack_p_brood_total)] <- 0
                
                #Now just calculate jack proportion of males
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                jack_p_brood_1_male <- matrix(nrow=n_brood, ncol=n_sim)
                #jacks sired by jacks
                jack_p_brood_2_male <- matrix(nrow=n_brood, ncol=n_sim)
                #total jacks
                jack_p_brood_total_male <- matrix(nrow=n_brood, ncol=n_sim)
                total_return_male <- matrix(nrow=n_brood, ncol=n_sim)
                for(y in 1:(n_brood)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_return_male[y,l] <- 0
                    for(a in 1:n_age){
                      for(j in 1:n_sire){
                        #conditional to avoid errors from empty lists
                        if(is.null(dim(mat_list[[y+a]][[a]][[2]][[j]][[l]]))==FALSE){
                          total_return_male[y,l] <- total_return_male[y,l] + nrow(mat_list[[y+a]][[a]][[2]][[j]][[l]]) 
                        }
                      }
                    }
                    jack_p_brood_1_male[y,l] <- nrow(mat_list[[y+3]][[3]][[2]][[1]][[l]])/total_return_male[y,l]
                    jack_p_brood_2_male[y,l] <- nrow(mat_list[[y+3]][[3]][[2]][[2]][[l]])/total_return_male[y,l]
                    jack_p_brood_total_male[y,l] <- (nrow(mat_list[[y+3]][[3]][[2]][[1]][[l]]) + nrow(mat_list[[y+3]][[3]][[2]][[2]][[l]]))/total_return_male[y,l]
                  }
                }
                
                #Plot time series of jack proportion across sire types and total, mean and 95% quantiles across simulations
                #total jack proportion
                plot(apply(jack_p_brood_total, 1, mean, na.rm=TRUE), type='l', col='black', ylim=c(0, max(apply(jack_p_brood_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                       'Jack proportion', main= 'Total jack prop (brood year)', cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_brood), rev(1:n_brood))
                yy <- c(apply(jack_p_brood_total, 1, cred)[1,], rev(apply(jack_p_brood_total, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by jacks
                points(apply(jack_p_brood_2, 1, mean, na.rm=TRUE), type='l', col='red', lwd=1.5)
                yy <- c(apply(jack_p_brood_2, 1, cred)[1,], rev(apply(jack_p_brood_2, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(100,1,1, max=250, alpha=90), border=FALSE)
                #jacks sired by non-jacks
                points(apply(jack_p_brood_1, 1, mean, na.rm=TRUE), type='l', col='blue', lwd=1.5)
                yy <- c(apply(jack_p_brood_1, 1, cred)[1,], rev(apply(jack_p_brood_1, 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(1,1,100, max=250, alpha=90), border=FALSE)
                legend(x=100, y=max(apply(jack_p_brood_total, 1, cred)[2,]), legend=c('jack-sired', 'hooknose-sired', 'total jacks'), cex=0.65, col=c('red', 'blue', 'black'), lty=1, bty='n')
                
                # #plot relationship between jacks in spawners and jacks in the return (time series)
                # plot(apply(jack_p_brood_total, 1, mean, na.rm=TRUE), type='l', col='red', ylim=c(0, max(apply(jack_p_spawn_total, 1, cred)[2,])), lwd=1.5, xlab='Year', ylab=
                #        'Jack proportion', main='Jack prop (brood and return)', cex.main=0.85, cex=0.85)
                # mtext(side=1, 'Year', cex=0.85, line=2)
                # mtext(side=2, 'Proportion', cex=0.85, line=2)
                # points(apply(jack_p_spawn_total, 1, mean, na.rm=TRUE), type='l', col='blue')
                # legend(x=100, y=max(apply(jack_p_spawn_total, 1, cred)[2,]), legend=c('jacks in spawners', 'jacks in return'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #plot time series of the population-wide age comp over time
                n_brood <- n_year-n_age
                #male brood year age comp
                age_comp_brood_male <- array(NA, dim=c(n_brood, n_sim, n_age))
                #female brood year age comp
                age_comp_brood_female <- array(NA, dim=c(n_brood, n_sim, n_age))
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                total_return_male <- matrix(nrow=n_brood, ncol=n_sim)
                total_return_female <- matrix(nrow=n_brood, ncol=n_sim)
                
                for(y in 1:(n_brood)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_return_male[y,l] <- 0
                    total_return_female[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(s==1){
                            if(is.null(dim(mat_list[[y+a]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_return_female[y,l] <- total_return_female[y,l] + nrow(mat_list[[y+a]][[a]][[s]][[j]][[l]]) 
                            }
                          }
                          if(s==2){
                            if(is.null(dim(mat_list[[y+a]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_return_male[y,l] <- total_return_male[y,l] + nrow(mat_list[[y+a]][[a]][[s]][[j]][[l]]) 
                            }
                          } 
                        }
                      }
                    }
                  }
                }  
                for(y in 1:(n_brood)){
                  for(l in 1:n_sim){
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        if(s==2){
                          if(is.null(dim(mat_list[[y+a]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_brood_male[y,l,a] <- (nrow(mat_list[[y+a]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y+a]][[a]][[s]][[2]][[l]]))/total_return_male[y,l]
                          }
                        }
                        if(s==1){
                          if(is.null(dim(mat_list[[y+a]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_brood_female[y,l,a] <- (nrow(mat_list[[y+a]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y+a]][[a]][[s]][[2]][[l]]))/total_return_female[y,l]
                          }
                        }
                      }
                    }
                  }
                }
                
                #replace NaNs in the age composition with 0s
                age_comp_brood_female[is.na(age_comp_brood_female)] <- 0
                age_comp_brood_male[is.na(age_comp_brood_male)] <- 0
                
                #plot time series (females)
                col.vec <- brewer.pal(n=(n_age), name='Blues')
                plot(apply(age_comp_brood_female[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='female age (brood year)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_brood), rev(1:n_brood))
                yy <- c(apply(age_comp_brood_female[,,3], 1, cred)[1,], rev(apply(age_comp_brood_female[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_brood_female[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='female age over time (brood year)')
                  xx <- c((1:n_brood), rev(1:n_brood))
                  yy <- c(apply(age_comp_brood_female[,,a], 1, cred)[1,], rev(apply(age_comp_brood_female[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #plot time series (males)
                plot(apply(age_comp_brood_male[,,3], 1, mean), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='male age (brood year)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_brood), rev(1:n_brood))
                yy <- c(apply(age_comp_brood_male[,,3], 1, cred)[1,], rev(apply(age_comp_brood_male[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_brood_male[,,a], 1, mean), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='male age over time (brood year)')
                  xx <- c((1:n_brood), rev(1:n_brood))
                  yy <- c(apply(age_comp_brood_male[,,a], 1, cred)[1,], rev(apply(age_comp_brood_male[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #age comp by esc year (escapement)
                n_esc <- n_year
                #male esc year age comp
                age_comp_esc_male <- array(NA, dim=c(n_esc, n_sim, n_age))
                #female esc year age comp
                age_comp_esc_female <- array(NA, dim=c(n_esc, n_sim, n_age))
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total spawning population in a given year
                total_spawner_male <- matrix(nrow=n_esc, ncol=n_sim)
                total_spawner_female <- matrix(nrow=n_esc, ncol=n_sim)
                
                for(y in 1:(n_esc-1)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_spawner_male[y,l] <- 0
                    total_spawner_female[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(s==1){
                            if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_spawner_female[y,l] <- total_spawner_female[y,l] + nrow(mat_list_esc[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          }
                          if(s==2){
                            if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_spawner_male[y,l] <- total_spawner_male[y,l] + nrow(mat_list_esc[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          } 
                        }
                      }
                    }
                  }
                }  
                
                for(y in 1:(n_esc-1)){
                  for(l in 1:n_sim){
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        if(s==2){
                          if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_esc_male[y,l,a] <- (nrow(mat_list_esc[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_esc[[y]][[a]][[s]][[2]][[l]]))/total_spawner_male[y,l]
                          }
                        }
                        if(s==1){
                          if(is.null(dim(mat_list_esc[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_esc_female[y,l,a] <- (nrow(mat_list_esc[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_esc[[y]][[a]][[s]][[2]][[l]]))/total_spawner_female[y,l]
                          }
                        }
                      }
                    }
                  }
                }
                
                #replace NaNs in the age composition with 0s
                age_comp_esc_female[is.na(age_comp_esc_female)] <- 0
                age_comp_esc_male[is.na(age_comp_esc_male)] <- 0
                
                
                #plot time series (females)
                col.vec <- brewer.pal(n=(n_age), name='Blues')
                plot(apply(age_comp_esc_female[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='female escapement age comp',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_esc), rev(1:n_esc))
                yy <- c(apply(age_comp_esc_female[,,3], 1, cred)[1,], rev(apply(age_comp_esc_female[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_esc_female[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='female age over time in escapement (return year)')
                  xx <- c((1:n_esc), rev(1:n_esc))
                  yy <- c(apply(age_comp_esc_female[,,a], 1, cred)[1,], rev(apply(age_comp_esc_female[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #plot time series (males)
                plot(apply(age_comp_esc_male[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='male escapement age comp',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_esc), rev(1:n_esc))
                yy <- c(apply(age_comp_esc_male[,,3], 1, cred)[1,], rev(apply(age_comp_esc_male[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_esc_male[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='male age over time in escapement (return year)')
                  xx <- c((1:n_esc), rev(1:n_esc))
                  yy <- c(apply(age_comp_esc_male[,,a], 1, cred)[1,], rev(apply(age_comp_esc_male[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #age comp by catch year (harvest)
                n_return <- n_year
                #male catch year age comp
                age_comp_catch_male <- array(NA, dim=c(n_return, n_sim, n_age))
                #female catch year age comp
                age_comp_catch_female <- array(NA, dim=c(n_return, n_sim, n_age))
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total catching population in a given year
                total_catch_male <- matrix(nrow=n_return, ncol=n_sim)
                total_catch_female <- matrix(nrow=n_return, ncol=n_sim)
                
                for(y in burnin:(n_return-1)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_catch_male[y,l] <- 0
                    total_catch_female[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(s==1){
                            if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_catch_female[y,l] <- total_catch_female[y,l] + nrow(mat_list_catch[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          }
                          if(s==2){
                            if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_catch_male[y,l] <- total_catch_male[y,l] + nrow(mat_list_catch[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          } 
                        }
                      }
                    }
                  }
                }  
                
                for(y in burnin:(n_return-1)){
                  for(l in 1:n_sim){
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        if(s==2){
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_catch_male[y,l,a] <- (nrow(mat_list_catch[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y]][[a]][[s]][[2]][[l]]))/total_catch_male[y,l]
                          }
                        }
                        if(s==1){
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_catch_female[y,l,a] <- (nrow(mat_list_catch[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y]][[a]][[s]][[2]][[l]]))/total_catch_female[y,l]
                          }
                        }
                      }
                    }
                  }
                }
                
                #replace NaNs in the age composition with 0s
                age_comp_catch_female[is.na(age_comp_catch_female)] <- 0
                age_comp_catch_male[is.na(age_comp_catch_male)] <- 0
                
                
                #plot time series (females)
                col.vec <- brewer.pal(n=(n_age), name='Blues')
                plot(apply(age_comp_catch_female[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='female harvest age comp',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_catch_female[,,3], 1, cred)[1,], rev(apply(age_comp_catch_female[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_catch_female[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='female age over time in harvest (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_catch_female[,,a], 1, cred)[1,], rev(apply(age_comp_catch_female[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #plot time series (males)
                plot(apply(age_comp_catch_male[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='male harvest age comp',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_catch_male[,,3], 1, cred)[1,], rev(apply(age_comp_catch_male[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_catch_male[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='male age over time in harvest (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_catch_male[,,a], 1, cred)[1,], rev(apply(age_comp_catch_male[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #Age-specific exploitation rates for males and females
                n_return <- n_year
                #male catch year age comp
                age_comp_catch_male <- array(NA, dim=c(n_return, n_sim, n_age))
                #female catch year age comp
                age_comp_catch_female <- array(NA, dim=c(n_return, n_sim, n_age))
                
                for(y in burnin:(n_return-1)){
                  for(l in 1:n_sim){
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        if(s==2){
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_catch_male[y,l,a] <- (nrow(mat_list_catch[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y]][[a]][[s]][[2]][[l]]))/(nrow(mat_list[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y]][[a]][[s]][[2]][[l]]))
                          }
                        }
                        if(s==1){
                          if(is.null(dim(mat_list_catch[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_catch_female[y,l,a] <- (nrow(mat_list_catch[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list_catch[[y]][[a]][[s]][[2]][[l]]))/(nrow(mat_list[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y]][[a]][[s]][[2]][[l]]))
                          }
                        }
                      }
                    }
                  }
                }
                
                #replace NaNs in the age composition with 0s
                age_comp_catch_female[is.na(age_comp_catch_female)] <- 0
                age_comp_catch_male[is.na(age_comp_catch_male)] <- 0
                
                
                #plot time series (females)
                col.vec <- brewer.pal(n=(n_age), name='Blues')
                plot(apply(age_comp_catch_female[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'exploitation rate', main='harvest rate at age (female)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Exploitation rate', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_catch_female[,,3], 1, cred)[1,], rev(apply(age_comp_catch_female[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_catch_female[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                           'exploitation rate', main='exploitation rates at age for females (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_catch_female[,,a], 1, cred)[1,], rev(apply(age_comp_catch_female[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #plot time series (males)
                plot(apply(age_comp_catch_male[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'exploitation rate', main='harvest rate at age (male)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Exploitation rate', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_catch_male[,,3], 1, cred)[1,], rev(apply(age_comp_catch_male[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_catch_male[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                           'exploitation rate', main='exploitation rates at age for males (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_catch_male[,,a], 1, cred)[1,], rev(apply(age_comp_catch_male[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #age comp by return year
                n_return <- n_year
                #male return year age comp
                age_comp_return_male <- array(NA, dim=c(n_return, n_sim, n_age))
                #female return year age comp
                age_comp_return_female <- array(NA, dim=c(n_return, n_sim, n_age))
                
                #loop through years and simulations --> produces a year x simulation # matrix of jacks abundance
                #matrix to store the total runing population in a given year
                total_run_male <- matrix(nrow=n_return, ncol=n_sim)
                total_run_female <- matrix(nrow=n_return, ncol=n_sim)
                
                for(y in 1:(n_return)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_run_male[y,l] <- 0
                    total_run_female[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(s==1){
                            if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_run_female[y,l] <- total_run_female[y,l] + nrow(mat_list[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          }
                          if(s==2){
                            if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_run_male[y,l] <- total_run_male[y,l] + nrow(mat_list[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          } 
                        }
                      }
                    }
                  }
                }  
                
                for(y in 1:(n_return)){
                  for(l in 1:n_sim){
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        if(s==2){
                          if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_return_male[y,l,a] <- (nrow(mat_list[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y]][[a]][[s]][[2]][[l]]))/total_run_male[y,l]
                          }
                        }
                        if(s==1){
                          if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                            age_comp_return_female[y,l,a] <- (nrow(mat_list[[y]][[a]][[s]][[1]][[l]]) + nrow(mat_list[[y]][[a]][[s]][[2]][[l]]))/total_run_female[y,l]
                          }
                        }
                      }
                    }
                  }
                }
                
                #replace NaNs in the age composition with 0s
                age_comp_return_female[is.na(age_comp_return_female)] <- 0
                age_comp_return_male[is.na(age_comp_return_male)] <- 0
                
                
                #plot time series (females)
                col.vec <- brewer.pal(n=(n_age), name='Blues')
                plot(apply(age_comp_return_female[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='total female age (return year)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_return_female[,,3], 1, cred)[1,], rev(apply(age_comp_return_female[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_return_female[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='female age over time (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_return_female[,,a], 1, cred)[1,], rev(apply(age_comp_return_female[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #plot time series (males)
                plot(apply(age_comp_return_male[,,3], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylim=c(0,1.25), lwd=1.5, xlab='Year', ylab=
                       'Proportion-at-age', main='total male age (return year)',cex.main=0.85, cex=0.85)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Proportion', cex=0.85, line=2)
                xx <- c((1:n_return), rev(1:n_return))
                yy <- c(apply(age_comp_return_male[,,3], 1, cred)[1,], rev(apply(age_comp_return_male[,,3], 1, cred)[2,]))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=200), border=FALSE)
                for(a in 4:n_age){
                  points(apply(age_comp_return_male[,,a], 1, mean, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), ylim=c(0,1), lwd=1.5, xlab='Year', ylab=
                           'Proportion-at-age', main='male age over time (return year)')
                  xx <- c((1:n_return), rev(1:n_return))
                  yy <- c(apply(age_comp_return_male[,,a], 1, cred)[1,], rev(apply(age_comp_return_male[,,a], 1, cred)[2,]))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=200), border=FALSE)
                }
                legend(x=100, y=1.35, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                
                #Plot y_0 valus of the PMRN over time
                for(y in 1:(n_return)){
                  for(l in 1:n_sim){
                    #initiate calculation of total spwaners at 0
                    total_spawner_male[y,l] <- 0
                    total_spawner_female[y,l] <- 0
                    for(a in 1:n_age){
                      for(s in 1:n_sex){
                        for(j in 1:n_sire){
                          #conditional to avoid errors from empty lists
                          if(s==1){
                            if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_spawner_female[y,l] <- total_spawner_female[y,l] + nrow(mat_list[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          }
                          if(s==2){
                            if(is.null(dim(mat_list[[y]][[a]][[s]][[j]][[l]]))==FALSE){
                              total_spawner_male[y,l] <- total_spawner_male[y,l] + nrow(mat_list[[y]][[a]][[s]][[j]][[l]]) 
                            }
                          } 
                        }
                      }
                    }
                  }
                }
                
                #calculate trait values given to recruits over time
                #y_0
                #median for males and females
                y_0_recruit_male_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                y_0_recruit_female_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                #lower 95% bounds for males and females
                y_0_recruit_male_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                y_0_recruit_female_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                #upper 95% bounds for males and females
                y_0_recruit_male_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                y_0_recruit_female_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                
                #loop through years, simulations, and sire types
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    for(j in 1:n_sire){
                      #calculate the median trait value of male recruits
                      y_0_recruit_male_mu[y,l,j] <- median(imm_list[[y]][[1]][[2]][[j]][[l]][,1], na.rm=TRUE)
                      #calculate the median trait value of female recruits
                      y_0_recruit_female_mu[y,l,j] <- median(imm_list[[y]][[1]][[1]][[j]][[l]][,1], na.rm=TRUE)
                      #calculate the lower bound of the trait value of male recruits
                      y_0_recruit_male_lower[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,1])[1]
                      #calculate the lower bound of the trait value of mfeale recruits
                      y_0_recruit_female_lower[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,1])[1]
                      #calculate the upper bound of the trait value of male recruits
                      y_0_recruit_male_upper[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,1])[2]
                      #calculate the upper bound of the trait value of female recruits
                      y_0_recruit_female_upper[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,1])[2]
                    }
                  }
                }
                
                #plot the median of the median trait value among simulations, as well as the median upper and lower bounds among simulations
                #males
                plot(apply(y_0_recruit_male_mu[,,1], 1, median, na.rm=TRUE), ylim=c(200,1800), type='l', col='blue', ylab='Recruit y_0 trait value', xlab='year', main='male PMRN intercept', cex.main=0.85, cex=0.85)
                points(apply(y_0_recruit_male_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(y_0_recruit_male_lower[,,1], 1, median, na.rm=TRUE), rev(apply(y_0_recruit_male_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(y_0_recruit_male_lower[,,2], 1, median, na.rm=TRUE), rev(apply(y_0_recruit_male_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'PMRN intercept', cex=0.85, line=2)
                legend(x=100, y=450, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #females
                plot(apply(y_0_recruit_female_mu[,,1], 1, median, na.rm=TRUE), ylim=c(200,1800), type='l', col='blue', ylab='Recruit y_0 trait value', xlab='year', main='female PMRN intercept', cex.main=0.85, cex=0.85)
                points(apply(y_0_recruit_female_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(y_0_recruit_female_lower[,,1], 1, median, na.rm=TRUE), rev(apply(y_0_recruit_female_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(y_0_recruit_female_lower[,,2], 1, median, na.rm=TRUE), rev(apply(y_0_recruit_female_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'PMRN intercept', cex=0.85, line=2)
                legend(x=100, y=450, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #Linfinity
                #median for males and females
                L_inf_recruit_male_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                L_inf_recruit_female_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                #lower 95% bounds for males and females
                L_inf_recruit_male_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                L_inf_recruit_female_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                #upper 95% bounds for males and females
                L_inf_recruit_male_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                L_inf_recruit_female_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                
                #loop through years, simulations, and sire types
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    for(j in 1:n_sire){
                      #calculate the median trait value of male recruits
                      L_inf_recruit_male_mu[y,l,j] <- median(imm_list[[y]][[1]][[2]][[j]][[l]][,2], na.rm=TRUE)
                      #calculate the median trait value of female recruits
                      L_inf_recruit_female_mu[y,l,j] <- median(imm_list[[y]][[1]][[1]][[j]][[l]][,2], na.rm=TRUE)
                      #calculate the lower bound of the trait value of male recruits
                      L_inf_recruit_male_lower[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,2])[1]
                      #calculate the lower bound of the trait value of mfeale recruits
                      L_inf_recruit_female_lower[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,2])[1]
                      #calculate the upper bound of the trait value of male recruits
                      L_inf_recruit_male_upper[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,2])[2]
                      #calculate the upper bound of the trait value of female recruits
                      L_inf_recruit_female_upper[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,2])[2]
                    }
                  }
                }
                
                #plot the median of the median trait value among simulations, as well as the median upper and lower bounds among simulations
                #males
                plot(apply(L_inf_recruit_male_mu[,,1], 1, median, na.rm=TRUE), ylim=c(200,1000), type='l', col='blue', ylab='Recruit L infinity trait value', xlab='year', main='male asymptotic length', cex.main=0.85, cex=0.85)
                points(apply(L_inf_recruit_male_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(L_inf_recruit_male_lower[,,1], 1, median, na.rm=TRUE), rev(apply(L_inf_recruit_male_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(L_inf_recruit_male_lower[,,2], 1, median, na.rm=TRUE), rev(apply(L_inf_recruit_male_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Asymptotic length', cex=0.85, line=2)
                legend(x=100, y=400, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #females
                plot(apply(L_inf_recruit_female_mu[,,1], 1, median, na.rm=TRUE), ylim=c(200,1000), type='l', col='blue', ylab='Recruit L infinity trait value', xlab='year', main='female asymptotic length', cex.main=0.85, cex=0.85)
                points(apply(L_inf_recruit_female_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(L_inf_recruit_female_lower[,,1], 1, median, na.rm=TRUE), rev(apply(L_inf_recruit_female_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(L_inf_recruit_female_lower[,,2], 1, median, na.rm=TRUE), rev(apply(L_inf_recruit_female_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Asymptotic length', cex=0.85, line=2)
                legend(x=100, y=400, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #Brody growth coefficient
                #median for males and females
                K_brody_recruit_male_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                K_brody_recruit_female_mu <- array(NA, dim=c(n_year, n_sim, n_sire))
                #lower 95% bounds for males and females
                K_brody_recruit_male_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                K_brody_recruit_female_lower <- array(NA, dim=c(n_year, n_sim, n_sire))
                #upper 95% bounds for males and females
                K_brody_recruit_male_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                K_brody_recruit_female_upper <- array(NA, dim=c(n_year, n_sim, n_sire))
                
                #loop through years, simulations, and sire types
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    for(j in 1:n_sire){
                      #calculate the median trait value of male recruits
                      K_brody_recruit_male_mu[y,l,j] <- median(imm_list[[y]][[1]][[2]][[j]][[l]][,3], na.rm=TRUE)
                      #calculate the median trait value of female recruits
                      K_brody_recruit_female_mu[y,l,j] <- median(imm_list[[y]][[1]][[1]][[j]][[l]][,3], na.rm=TRUE)
                      #calculate the lower bound of the trait value of male recruits
                      K_brody_recruit_male_lower[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,3])[1]
                      #calculate the lower bound of the trait value of mfeale recruits
                      K_brody_recruit_female_lower[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,3])[1]
                      #calculate the upper bound of the trait value of male recruits
                      K_brody_recruit_male_upper[y,l,j] <- cred(imm_list[[y]][[1]][[2]][[j]][[l]][,3])[2]
                      #calculate the upper bound of the trait value of female recruits
                      K_brody_recruit_female_upper[y,l,j] <- cred(imm_list[[y]][[1]][[1]][[j]][[l]][,3])[2]
                    }
                  }
                }
                
                #plot the median of the median trait value among simulations, as well as the median upper and lower bounds among simulations
                #males
                plot(apply(K_brody_recruit_male_mu[,,1], 1, median, na.rm=TRUE), ylim=c(0,1), type='l', col='blue', ylab='Recruit K trait value', xlab='year', main='male growth rate', cex.main=0.85, cex=0.85)
                points(apply(K_brody_recruit_male_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(K_brody_recruit_male_lower[,,1], 1, median, na.rm=TRUE), rev(apply(K_brody_recruit_male_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(K_brody_recruit_male_lower[,,2], 1, median, na.rm=TRUE), rev(apply(K_brody_recruit_male_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Growth rate', cex=0.85, line=2)
                legend(x=100, y=0.9, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
                
                #females
                plot(apply(K_brody_recruit_female_mu[,,1], 1, median, na.rm=TRUE), ylim=c(0,1), type='l', col='blue', ylab='Recruit K trait value', xlab='year', main='female growth rate', cex.main=0.85, cex=0.85)
                points(apply(K_brody_recruit_female_mu[,,2], 1, median, na.rm=TRUE), type='l', col='red')
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(K_brody_recruit_female_lower[,,1], 1, median, na.rm=TRUE), rev(apply(K_brody_recruit_female_upper[,,1], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='blue')[1], col2rgb(col='blue')[2], col2rgb(col='blue')[3], max=255, alpha=200), border=FALSE)
                yy <- c(apply(K_brody_recruit_female_lower[,,2], 1, median, na.rm=TRUE), rev(apply(K_brody_recruit_female_upper[,,2], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col='red')[1], col2rgb(col='red')[2], col2rgb(col='red')[3], max=255, alpha=100), border=FALSE)
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Growth rate', cex=0.85, line=2)
                legend(x=100, y=0.9, legend=c('jack-sired', 'hooknose-sired'), cex=0.65, col=c('red', 'blue'), lty=1, bty='n')
               
                #Length-at-age (spawners)
                #median for males and females
                length_at_age_recruit_male_mu <- array(NA, dim=c(n_year, n_sim, n_age))
                length_at_age_recruit_female_mu <- array(NA, dim=c(n_year, n_sim, n_age))
                #lower 95% bounds for males and females
                length_at_age_recruit_male_lower <- array(NA, dim=c(n_year, n_sim, n_age))
                length_at_age_recruit_female_lower <- array(NA, dim=c(n_year, n_sim, n_age))
                #upper 95% bounds for males and females
                length_at_age_recruit_male_upper <- array(NA, dim=c(n_year, n_sim, n_age))
                length_at_age_recruit_female_upper <- array(NA, dim=c(n_year, n_sim, n_age))
                
                #loop through years, simulations, and ages types
                for(y in 1:n_year){
                  for(l in 1:n_sim){
                    for(a in 3:n_age){
                      #calculate the median trait value of male spawners
                      length_at_age_recruit_male_mu[y,l,a] <- median(mat_list[[y]][[a]][[2]][[1]][[l]][,4], mat_list[[y]][[a]][[2]][[2]][[l]][,4], na.rm=TRUE)
                      #calculate the median trait value of female spawners
                      length_at_age_recruit_female_mu[y,l,a] <- median(mat_list[[y]][[a]][[1]][[1]][[l]][,4], mat_list[[y]][[a]][[1]][[2]][[l]][,4], na.rm=TRUE)
                      #calculate the lower bound of the trait value of male spawners
                      length_at_age_recruit_male_lower[y,l,a] <- cred(c(mat_list[[y]][[a]][[2]][[1]][[l]][,4], mat_list[[y]][[a]][[2]][[2]][[l]][,4]))[1]
                      #calculate the lower bound of the trait value of female spawners
                      length_at_age_recruit_female_lower[y,l,a] <- cred(c(mat_list[[y]][[a]][[1]][[1]][[l]][,4], mat_list[[y]][[a]][[1]][[2]][[l]][,4]))[1]
                      #calculate the upper bound of the trait value of male spawners
                      length_at_age_recruit_male_upper[y,l,a] <- cred(c(mat_list[[y]][[a]][[2]][[1]][[l]][,4], mat_list[[y]][[a]][[2]][[2]][[l]][,4]))[2]
                      #calculate the upper bound of the trait value of female spawners
                      length_at_age_recruit_female_upper[y,l,a] <- cred(c(mat_list[[y]][[a]][[1]][[1]][[l]][,4], mat_list[[y]][[a]][[1]][[2]][[l]][,4]))[2]
                    }
                  }
                }
                
                #plot the median of the median trait value among simulations, as well as the median upper and lower bounds among simulations
                #males
                plot(apply(length_at_age_recruit_male_mu[,,3], 1, median, na.rm=TRUE), ylim=c(0,700), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylab='Length-at-age (spawners)', xlab='year', main='males', cex.main=0.85)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(length_at_age_recruit_male_lower[,,3], 1, median, na.rm=TRUE), rev(apply(length_at_age_recruit_male_upper[,,3], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), border=FALSE)
                for (a in 4:n_age){
                  points(apply(length_at_age_recruit_male_mu[,,a], 1, median, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255))
                  xx <- c((1:n_year), rev(1:n_year))
                  yy <- c(apply(length_at_age_recruit_male_lower[,,a], 1, median, na.rm=TRUE), rev(apply(length_at_age_recruit_male_upper[,,a], 1, median, na.rm=TRUE)))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), border=FALSE)
                }
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Size-at-age', cex=0.85, line=2)
                legend(x=100, y=200, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                #females
                plot(apply(length_at_age_recruit_female_mu[,,3], 1, median, na.rm=TRUE), ylim=c(0,700), type='l', col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), ylab='Length-at-age (spawners)', xlab='year', main='females', cex.main=0.85)
                xx <- c((1:n_year), rev(1:n_year))
                yy <- c(apply(length_at_age_recruit_female_lower[,,3], 1, median, na.rm=TRUE), rev(apply(length_at_age_recruit_female_upper[,,3], 1, median, na.rm=TRUE)))
                polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[3])[1], col2rgb(col=col.vec[3])[2], col2rgb(col=col.vec[3])[3], max=255, alpha=255), border=FALSE)
                for (a in 4:n_age){
                  points(apply(length_at_age_recruit_female_mu[,,a], 1, median, na.rm=TRUE), type='l', col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255))
                  xx <- c((1:n_year), rev(1:n_year))
                  yy <- c(apply(length_at_age_recruit_female_lower[,,a], 1, median, na.rm=TRUE), rev(apply(length_at_age_recruit_female_upper[,,a], 1, median, na.rm=TRUE)))
                  polygon(x=xx, y=yy, col=rgb(col2rgb(col=col.vec[a])[1], col2rgb(col=col.vec[a])[2], col2rgb(col=col.vec[a])[3], max=255, alpha=255), border=FALSE)
                }
                mtext(side=1, 'Year', cex=0.85, line=2)
                mtext(side=2, 'Size-at-age', cex=0.85, line=2)
                legend(x=100, y=200, legend=c('Ocean age 1', 'Ocean age 2', 'Ocean age 3'), cex=0.65, col=col.vec[c(3:n_age)], lty=1, bty='n')
                dev.off()
              }
              jack_mat[[r]][[w]][x,h]   <- median(apply(jack_p_brood_total[(n_year-250):(n_year-n_age-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
              jack_mat_return[[r]][[w]][x,h]   <- median(apply(jack_p_run_total[(n_year-250):(n_year-n_age-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
              jack_mat_male[[r]][[w]][x,h]  <- median(apply(jack_p_brood_total_male[(n_year-250):(n_year-n_age-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
              jack_mat_return_male[[r]][[w]][x,h]   <- median(apply(jack_p_run_total_male[(n_year-250):(n_year-n_age-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
              catch_mat[[r]][[w]][x,h]   <-  median(apply(total_catch[(n_year-250):(n_year-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
              u_real_mat[[r]][[w]][x,h]   <-  median(apply(u_real[502:511,],2, mean, na.rm=TRUE),na.rm=TRUE)
              non_jack_catch_mat[[r]][[w]][x,h]   <-  median(apply(non_jack_catch[(n_year-250):(n_year-1),],2, mean, na.rm=TRUE),na.rm=TRUE)
            }
          }
        }
      }
    }
  }

saveRDS(jack_mat, file="jack_mat_h.RDS")
saveRDS(jack_mat_return, file="jack_mat_return_h.RDS")
saveRDS(jack_mat_male, file="jack_mat_male_h.RDS")
saveRDS(jack_mat_return_male, file="jack_mat_return_male_h.RDS")
saveRDS(catch_mat, file="catch_mat_h.RDS")
saveRDS(u_real_mat, file="u_real_mat_h.RDS")
saveRDS(non_jack_catch_mat, file="non_jack_catch_mat_h.RDS")

