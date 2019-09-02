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
library(ggplot2)

###
# Parallelize Stan.
###

no_cores <- min(3, (parallel::detectCores()-1))
options(mc.cores = (no_cores))

###
# Starting parameters.
###

n_chains <- 3
n_iter <- 2000
n_warmup <- 1000

zeta <- -3.5
tau <- c(3.8, 2.9)

phis <- c(1, 2, 3)
n_phis <- length(phis)
psi_CMs <- c(1, 2, 3); psi_CPs <- c(1, 2, 3)
n_psis <- length(psi_CMs)

mads <- list(list(list(list())))
n_datasets <- 3
for(i in 1:n_datasets){
  for(j in 1:n_phis){
    for(k in 1:n_psis){
      
      ###
      # Generate simulated data.
      ###
      
      N <- 1000
      
      ME <- rbinom(N, 1, 0.1)
      PE <- rbinom(N, 1, 0.1)
      X <- cbind(ME, PE)
      K <- ncol(X)
      
      gamma <- (zeta + X %*% tau)
      phi <- phis[j]
      trans_theta <- rnorm(N, gamma, phi)
      theta <- exp(trans_theta) / (1+exp(trans_theta))
      
      psi_CM <- psi_CMs[k]
      psi_CP <- psi_CPs[k]
      trans_CP <- rnorm(N, trans_theta, psi_CP)
      trans_CM <- rnorm(N, trans_theta, psi_CM)
      CP <- exp(trans_CP) / (1+exp(trans_CP))
      CM <- exp(trans_CM) / (1+exp(trans_CM))
      Y <- cbind(CP, CM)
      M <- ncol(Y)
      
      proxy_structure <- data.frame("regulator"="Regulator", "target_candidate"=paste0("TargetCandidate",1:N), "ME"=ME, "PE"=PE, "CM"=CM, "CP"=CP)
      prepared_data <- prepare_data(proxy_structure)
      
      ###
      # Full BINDER performance (MADs).
      ###
      
      model_binder <- run_binder(prepared_data, mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_phi=0, sigma_phi=3, mu_psi=c(0,0), sigma_psi=c(3,3), iter=n_iter, warmup=n_warmup, chains=n_chains, seed=1)
      theta_samples_binder <- rstan::extract(model_binder$model_object, "theta")[[1]]
      
      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_binder <- rep(0,N)
      for(n in 1:N){
        theta_sample <- as.numeric(theta_samples_binder[, n])
        theta_mad_binder <- mean(abs(theta[n] - theta_sample))
        theta_mads_binder[n] <- theta_mad_binder
      }
      mads[[as.character(i)]][[as.character(j)]][[as.character(k)]][["binder"]] <- theta_mads_binder
      
      ###
      # Non-auxiliary BINDER performance (MADs).
      ###
      
      model_non_auxiliary <- run_non_auxiliary_binder(prepared_data, mu_psi=c(0,0), sigma_psi=c(3,3), iter=n_iter, warmup=n_warmup, chains=n_chains, seed=1)
      theta_samples_non_auxiliary <- rstan::extract(model_non_auxiliary$model_object, "theta")[[1]]
      
      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_non_auxiliary <- rep(0,N)
      for(n in 1:N){
        theta_sample <- as.numeric(theta_samples_non_auxiliary[, n])
        theta_mad_non_auxiliary <- mean(abs(theta[n] - theta_sample))
        theta_mads_non_auxiliary[n] <- theta_mad_non_auxiliary
      }
      mads[[as.character(i)]][[as.character(j)]][[as.character(k)]][["non_auxiliary"]] <- theta_mads_non_auxiliary
      
      ###
      # Deterministic BINDER performance (MADs).
      ###
      
      model_deterministic <- run_deterministic_binder(prepared_data, mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_psi=c(0,0), sigma_psi=c(3,3), iter=n_iter, warmup=n_warmup, chains=n_chains, seed=1)
      theta_samples_deterministic <- rstan::extract(model_deterministic$model_object, "theta")[[1]]
      
      n_samples <- n_chains*(n_iter-n_warmup)
      theta_mads_deterministic <- rep(0,N)
      for(n in 1:N){
        theta_sample <- as.numeric(theta_samples_deterministic[, n])
        theta_mad_deterministic <- mean(abs(theta[n] - theta_sample))
        theta_mads_deterministic[n] <- theta_mad_deterministic
      }
      mads[[as.character(i)]][[as.character(j)]][[as.character(k)]][["deterministic"]] <- theta_mads_deterministic
      
    }
  }
}

df <- data.frame(model=factor(rep(c("BINDER", "Non-auxiliary", "Deterministic"), each=9), levels=c("Deterministic", "Non-auxiliary", "BINDER")), phi=factor(rep(rep(c("low", "mid", "high"), each=3), 3), levels=c("high", "mid", "low")), psi=factor(rep(rep(c("psi[low]", "psi[mid]", "psi[high]"), 3), 3), levels=c("psi[low]", "psi[mid]", "psi[high]")))
mean_mads <- sapply(1:3, function(z){ sapply(c("binder", "non_auxiliary", "deterministic"), function(w){ sapply(1:3, function(y){ sapply(1:3, function(x){ mean(mads[[as.character(z)]][[as.character(y)]][[as.character(x)]][[w]]) }) }) }) })

df$mean <- apply(mean_mads, 1, mean)
df$se <- apply(mean_mads, 1, sd)
df$se_lower <- df$mean + df$se
df$se_upper <- df$mean - df$se

base_palette <- rev(c("#2E0219","#CE8964","#EFC69B"))

ggplot(df) +
  geom_bar(aes(x=phi, y=mean, fill=model, group=factor(1:(n_phis*n_psis*3))), stat="identity", position=position_dodge(1), color="white", alpha=0.6) +
  geom_errorbar(aes(x=phi, ymin=se_lower, ymax=se_upper, color=model, group=factor(1:(n_phis*n_psis*3))), width=0.5, position=position_dodge(1), size=1) +
  xlab(expression(paste(phi))) +
  ylab(expression(paste("MAD"[theta]))) +
  facet_grid(psi ~ ., labeller=label_parsed) +
  theme_minimal() +
  theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=12), legend.text=element_text(size=10), strip.text=element_text(size=12)) +
  scale_fill_manual(name="Model", values=base_palette) +
  scale_color_manual(name="Model", values=base_palette) +
  coord_flip()

mean(df$mean[df$model == "BINDER"])
sd(df$mean[df$model == "BINDER"])

mean(df$mean[df$model == "Non-auxiliary"])
sd(df$mean[df$model == "Non-auxiliary"])

mean(df$mean[df$model == "Deterministic"])
sd(df$mean[df$model == "Deterministic"])