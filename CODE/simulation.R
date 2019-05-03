###
# Set seed for replication.
###

set.seed(1)

#########################
#########################

###
# Import required libraries.
###

library(BINDER)

###
# Parallelize Stan.
###

no_cores <- min(3,(parallel::detectCores()-1))
options(mc.cores = (no_cores))

###
# Choose parameters.
###

phis <- c(0.1,1.5,3)
n_phis <- length(phis)
psi_CPs <- c(0.1,1.5,3)
n_psi_CPs <- length(psi_CPs)
psi_CMs <- c(0.1,1.5,3)
n_psi_CMs <- length(psi_CMs)


###
# Generate simulated data.
###

sim_data <- list()

n_datasets <- 3
for(dataset in 1:n_datasets){
  sim_data[[as.character(dataset)]] <- list()
  iter <- 0
  for(first in 1:n_phis){
    for(second in 1:n_psi_CPs){

      iter <- iter+1

      sim_data[[as.character(dataset)]][[as.character(iter)]]$N <- 1000

      sim_data[[as.character(dataset)]][[as.character(iter)]]$ME <- rbinom(sim_data[[as.character(dataset)]][[as.character(iter)]]$N,1,0.01)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$PE <- rbinom(sim_data[[as.character(dataset)]][[as.character(iter)]]$N,1,0.01)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$X <- cbind(sim_data[[as.character(dataset)]][[as.character(iter)]]$ME, sim_data[[as.character(dataset)]][[as.character(iter)]]$PE)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$K <- dim(sim_data[[as.character(dataset)]][[as.character(iter)]]$X)[2]

      sim_data[[as.character(dataset)]][[as.character(iter)]]$zeta <- -1.1
      sim_data[[as.character(dataset)]][[as.character(iter)]]$tau <- c(1.8,0.9)

      sim_data[[as.character(dataset)]][[as.character(iter)]]$gamma <- sim_data[[as.character(dataset)]][[as.character(iter)]]$zeta + (sim_data[[as.character(dataset)]][[as.character(iter)]]$X %*% sim_data[[as.character(dataset)]][[as.character(iter)]]$tau)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$phi <- phis[first]
      sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_theta <- rnorm(sim_data[[as.character(dataset)]][[as.character(iter)]]$N, sim_data[[as.character(dataset)]][[as.character(iter)]]$gamma, sim_data[[as.character(dataset)]][[as.character(iter)]]$phi)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$theta <- exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_theta) / (1+exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_theta))

      sim_data[[as.character(dataset)]][[as.character(iter)]]$psi_CP <- psi_CPs[second]
      sim_data[[as.character(dataset)]][[as.character(iter)]]$psi_CM <- psi_CMs[second]

      sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CP <- rnorm(sim_data[[as.character(dataset)]][[as.character(iter)]]$N, sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_theta, sim_data[[as.character(dataset)]][[as.character(iter)]]$psi_CP)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CM <- rnorm(sim_data[[as.character(dataset)]][[as.character(iter)]]$N, sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_theta, sim_data[[as.character(dataset)]][[as.character(iter)]]$psi_CM)

      sim_data[[as.character(dataset)]][[as.character(iter)]]$CP <- exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CP) / (1+exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CP))
      sim_data[[as.character(dataset)]][[as.character(iter)]]$CM <- exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CM) / (1+exp(sim_data[[as.character(dataset)]][[as.character(iter)]]$trans_CM))

      sim_data[[as.character(dataset)]][[as.character(iter)]]$Y <- cbind(sim_data[[as.character(dataset)]][[as.character(iter)]]$CP, sim_data[[as.character(dataset)]][[as.character(iter)]]$CM)
      sim_data[[as.character(dataset)]][[as.character(iter)]]$M <- dim(sim_data[[as.character(dataset)]][[as.character(iter)]]$Y)[2]

      sim_data[[as.character(dataset)]][[as.character(iter)]]$proxy_structure <- data.frame("regulator"="Regulator", "target_candidate"=paste0("TargetCandidate",1:sim_data[[as.character(dataset)]][[as.character(iter)]]$N), "ME"=sim_data[[as.character(dataset)]][[as.character(iter)]]$ME, "PE"=sim_data[[as.character(dataset)]][[as.character(iter)]]$PE, "CM"=sim_data[[as.character(dataset)]][[as.character(iter)]]$CM, "CP"=sim_data[[as.character(dataset)]][[as.character(iter)]]$CP)

    }
  }
}

###
# Full BINDER performance (MADs).
###

all_theta_mads_binder <- c()

for(dataset in 1:n_datasets){
  iter <- 0
  for(first in 1:n_phis){
    for(second in 1:n_psi_CPs){
      
      iter <- iter+1
      
      n_chains <- 3
      n_iter <- 150
      n_warmup <- 100
      seed <- 1
      
      prepared_data <- prepare_data(sim_data[[as.character(dataset)]][[as.character(iter)]]$proxy_structure)
      model <- run_binder(prepared_data, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
      
      theta_samples_binder <- rstan::extract(model$model_object,"theta")[[1]]
      
      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_binder <- rep(0,sim_data[[as.character(dataset)]][[as.character(iter)]]$N) #rep(0,n_samples)
      for(n in 1:sim_data[[as.character(dataset)]][[as.character(iter)]]$N){ #for(n in 1:n_samples){
        theta_sample <- as.numeric(theta_samples_binder[, n]) #as.numeric(theta_samples_binder[n,])
        theta_mad_binder <- mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta[n] - theta_sample)) #mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta - theta_sample))
        theta_mads_binder[n] <- theta_mad_binder
      }
      sim_data[[as.character(dataset)]][[as.character(iter)]]$theta_mads_binder <- theta_mads_binder
      
      all_theta_mads_binder <- c(all_theta_mads_binder, theta_mads_binder)
      
    }
  }
}
  
dat <- rbind(sapply(1:n_datasets, function(y){ sapply(1:(n_phis*n_psi_CPs), function(x){ mean(sim_data[[as.character(y)]][[as.character(x)]]$theta_mads_binder) }) }))
dat <- data.frame(dat)
dat_binder <- t(cbind(sapply(1:nrow(dat), function(x){ current_mean <- mean(as.numeric(dat[x, ])); current_sd <- sd(as.numeric(dat[x, ])); cbind( current_mean, (current_mean-current_sd), (current_mean+current_sd) )  })))
dat_binder <- data.frame("Mean"=dat_binder[,1], "MeanMinusSD"=dat_binder[,2], "MeanPlusSD"=dat_binder[,3])

###
# Non-auxiliary BINDER performance (MADs).
###

all_theta_mads_non_auxiliary <- c()

for(dataset in 1:n_datasets){
  iter <- 0
  for(first in 1:n_phis){
    for(second in 1:n_psi_CPs){
      
      iter <- iter+1
      
      n_chains <- 3
      n_iter <- 150
      n_warmup <- 100
      seed <- 1
      
      prepared_data <- prepare_data(sim_data[[as.character(dataset)]][[as.character(iter)]]$proxy_structure)
      model <- run_non_auxiliary_binder(prepared_data, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)

      theta_samples_non_auxiliary <- rstan::extract(model$model_object,"theta")[[1]]

      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_non_auxiliary <- rep(0,sim_data[[as.character(dataset)]][[as.character(iter)]]$N) #rep(0,n_samples)
      for(n in 1:sim_data[[as.character(dataset)]][[as.character(iter)]]$N){ #for(n in 1:n_samples){
        theta_sample <- as.numeric(theta_samples_non_auxiliary[, n]) #as.numeric(theta_samples_non_auxiliary[n,])
        theta_mad_non_auxiliary <- mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta[n] - theta_sample)) #mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta - theta_sample))
        theta_mads_non_auxiliary[n] <- theta_mad_non_auxiliary
      }
      sim_data[[as.character(dataset)]][[as.character(iter)]]$theta_mads_non_auxiliary <- theta_mads_non_auxiliary
      
      all_theta_mads_non_auxiliary <- c(all_theta_mads_non_auxiliary, theta_mads_non_auxiliary)
      
    }
  }
}

dat <- rbind(sapply(1:n_datasets, function(y){ sapply(1:(n_phis*n_psi_CPs), function(x){ mean(sim_data[[as.character(y)]][[as.character(x)]]$theta_mads_non_auxiliary) }) }))
dat <- data.frame(dat)
dat_non_auxiliary <- t(cbind(sapply(1:nrow(dat), function(x){ current_mean <- mean(as.numeric(dat[x, ])); current_sd <- sd(as.numeric(dat[x, ])); cbind( current_mean, (current_mean-current_sd), (current_mean+current_sd) )  })))
dat_non_auxiliary <- data.frame("Mean"=dat_non_auxiliary[,1], "MeanMinusSD"=dat_non_auxiliary[,2] ,"MeanPlusSD"=dat_non_auxiliary[,3])


###
# deterministic BINDER performance (MADs).
###

all_theta_mads_deterministic <- c()

for(dataset in 1:n_datasets){
  iter <- 0
  for(first in 1:n_phis){
    for(second in 1:n_psi_CPs){
      
      iter <- iter+1
      
      n_chains <- 3
      n_iter <- 150
      n_warmup <- 100
      seed <- 1
      
      prepared_data <- prepare_data(sim_data[[as.character(dataset)]][[as.character(iter)]]$proxy_structure)
      model <- run_deterministic_binder(prepared_data, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
      
      theta_samples_deterministic <- rstan::extract(model$model_object,"theta")[[1]]
      
      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_deterministic <- rep(0,sim_data[[as.character(dataset)]][[as.character(iter)]]$N) #rep(0,n_samples)
      for(n in 1:sim_data[[as.character(dataset)]][[as.character(iter)]]$N){ #for(n in 1:n_samples){
        theta_sample <- as.numeric(theta_samples_deterministic[, n]) #as.numeric(theta_samples_deterministic[n,])
        theta_mad_deterministic <- mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta[n] - theta_sample)) #mean(abs(sim_data[[as.character(dataset)]][[as.character(iter)]]$theta - theta_sample))
        theta_mads_deterministic[n] <- theta_mad_deterministic
      }
      sim_data[[as.character(dataset)]][[as.character(iter)]]$theta_mads_deterministic <- theta_mads_deterministic
      
      all_theta_mads_deterministic <- c(all_theta_mads_deterministic, theta_mads_deterministic)
      
    }
  }
}

dat <- rbind(sapply(1:n_datasets, function(y){ sapply(1:(n_phis*n_psi_CPs), function(x){ mean(sim_data[[as.character(y)]][[as.character(x)]]$theta_mads_deterministic) }) }))
dat <- data.frame(dat)
dat_deterministic <- t(cbind(sapply(1:nrow(dat), function(x){ current_mean <- mean(as.numeric(dat[x, ])); current_sd <- sd(as.numeric(dat[x, ])); cbind( current_mean, (current_mean-current_sd), (current_mean+current_sd) )  })))
dat_deterministic <- data.frame("Mean"=dat_deterministic[,1], "MeanMinusSD"=dat_deterministic[,2] ,"MeanPlusSD"=dat_deterministic[,3])



###
# Generate results.
###

all_dat <- data.frame(rbind(as.matrix(dat_binder), as.matrix(dat_non_auxiliary), as.matrix(dat_deterministic)))

model_names <- c("BINDER","non_auxiliary","deterministic")
n_model_names <- length(model_names)

phi_names <- as.character(phis)
n_phi_names <- length(phi_names)

psi_names <- as.character(psi_CPs)
n_psi_names <- length(psi_names)

phi_g <- rep(rep(phi_names, each=n_psi_names), n_model_names)
psi_MP_g <- rep(rep(psi_names, n_phi_names), n_model_names)
model_g <- rep(c("BINDER","non-auxiliary","deterministic"), each=(n_phi_names*n_psi_names))

all_dat$phi_g <- factor(phi_g)
all_dat$psi_MP_g <- factor(psi_MP_g)
all_dat$model_g <- factor(model_g)
levels(all_dat$phi_g) <- c("low","mid","high")
all_dat$phi_g <- factor(all_dat$phi_g, levels=c("high","mid","low"))
levels(all_dat$psi_MP_g) <- c("psi[low]","psi[mid]","psi[high]")
all_dat$model_g <- factor(all_dat$model_g, levels=c("deterministic","non-auxiliary","BINDER"))

base_palette <- rev(c("#2E0219","#CE8964","#EFC69B"))

ggplot(all_dat) +
  geom_bar(aes(x=phi_g, y=all_dat$Mean, fill=model_g, group=factor(1:(n_model_names*n_phi_names*n_psi_names))), stat="identity", position=position_dodge(1), color="white", alpha=0.6) +
  geom_errorbar(aes(x=phi_g, ymin=all_dat$MeanMinusSD,ymax=all_dat$MeanPlusSD, color=model_g, group=factor(1:(n_model_names*n_phi_names*n_psi_names))), width=0.5, position=position_dodge(1), size=1) +
  xlab(expression(paste(phi))) +
  ylab(expression(paste("MAD"[theta]))) +
  scale_fill_manual(name="Model", values=base_palette) +
  scale_color_manual(name="Model", values=base_palette) +
  theme_minimal() +
  coord_flip() +
  facet_grid(psi_MP_g ~ ., labeller=label_parsed) +  theme(axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), strip.text=element_text(size=20))

mean(all_dat$Mean[all_dat$model_g=="BINDER"])
sd(all_dat$Mean[all_dat$model_g=="BINDER"])

mean(all_dat$Mean[all_dat$model_g=="non-auxiliary"])
sd(all_dat$Mean[all_dat$model_g=="non-auxiliary"])

mean(all_dat$Mean[all_dat$model_g=="deterministic"])
sd(all_dat$Mean[all_dat$model_g=="deterministic"])