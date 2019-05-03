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

#########################
#########################

###
# Parallelize Stan.
###

no_cores <- min(3, (parallel::detectCores()-1))
options(mc.cores=(no_cores))

#########################
#########################

###
# DATA.
###

m_abscessus_expression <- read.delim("DATA/m_abscessus/expression.txt", header=TRUE, row.names=1, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

m_abscessus_operons <- read.delim("DATA/m_abscessus/operons.opr", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
m_abscessus_operons <- m_abscessus_operons[order(m_abscessus_operons$Synonym), ]

keep <- intersect(m_abscessus_operons$Synonym, rownames(m_abscessus_expression))
m_abscessus_expression <- m_abscessus_expression[keep, ]
m_abscessus_operons <- m_abscessus_operons[m_abscessus_operons$Synonym %in% keep, ]
m_abscessus_coexpression <- compute_coexpression(m_abscessus_expression, keep)

all_proxy_regulons <- read.delim("DATA/m_abscessus/m_tuberculosis/ortholog_ME_PE.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
all_regulators <- unique(all_proxy_regulons$regulator)
n_regulators <- length(all_regulators)

#########################
#########################

###
# Run BINDER.
###

n_chains <- 3
n_iter <- 150
n_warmup <- 100
seed <- 1

#vague:
regulons_vague <- list()
for(n in 1:n_regulators){
  print(n)
  regulator <- all_regulators[n]
  print(regulator)
  proxy_regulon <- all_proxy_regulons[all_proxy_regulons$regulator == regulator, ]
  regulons_vague[[regulator]] <- binder(proxy_regulon, m_abscessus_coexpression, m_abscessus_operons, is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
}

final_regulons_vague <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_vague[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_vague[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_vague[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0){
    print(n)
    final_regulons_vague <- rbind( final_regulons_vague,
    data.frame(
     regulons_vague[[n]]$regulator,
     regulons_vague[[n]]$target_candidate,
     as.numeric(regulons_vague[[n]]$ME),
     as.numeric(regulons_vague[[n]]$PE),
     as.numeric(regulons_vague[[n]]$CM),
     as.numeric(regulons_vague[[n]]$CP),
     as.numeric(regulons_vague[[n]]$mean_theta),
     as.numeric(regulons_vague[[n]]$theta_interval[,1]),
     as.numeric(regulons_vague[[n]]$theta_interval[,2]),
     as.numeric(regulons_vague[[n]]$theta_interval[,3]),
     as.numeric(regulons_vague[[n]]$mean_gamma),
     as.numeric(regulons_vague[[n]]$gamma_interval[,1]),
     as.numeric(regulons_vague[[n]]$gamma_interval[,2]),
     as.numeric(regulons_vague[[n]]$mean_raw_gamma),
     as.numeric(regulons_vague[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_vague[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_vague[[n]]$mean_phi[1]),
     as.numeric(regulons_vague[[n]]$phi_interval[1]),
     as.numeric(regulons_vague[[n]]$phi_interval[2]),
     as.numeric(regulons_vague[[n]]$mean_zeta[1]),
     as.numeric(regulons_vague[[n]]$zeta_interval[1]),
     as.numeric(regulons_vague[[n]]$zeta_interval[2]),
     as.numeric(regulons_vague[[n]]$mean_tau[1]),
     as.numeric(regulons_vague[[n]]$tau_interval[1,1]),
     as.numeric(regulons_vague[[n]]$tau_interval[1,2]),
     as.numeric(regulons_vague[[n]]$mean_tau[2]),
     as.numeric(regulons_vague[[n]]$tau_interval[2,1]),
     as.numeric(regulons_vague[[n]]$tau_interval[2,2]),
     as.numeric(regulons_vague[[n]]$mean_psi[1]),
     as.numeric(regulons_vague[[n]]$psi_interval[1,1]),
     as.numeric(regulons_vague[[n]]$psi_interval[1,2]),
     as.numeric(regulons_vague[[n]]$mean_psi[2]),
     as.numeric(regulons_vague[[n]]$psi_interval[2,1]),
     as.numeric(regulons_vague[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_vague) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_vague, file="RESULTS/regulons_results/final_regulons.vague.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

##################################################
##################################################

###
# Run BINDER.
###

n_chains <- 3
n_iter <- 150
n_warmup <- 100
seed <- 1

#informed:
regulons_informed <- list()
for(n in 1:n_regulators){
  print(n)
  regulator <- all_regulators[n]
  print(regulator)
  proxy_regulon <- all_proxy_regulons[all_proxy_regulons$regulator == regulator, ]
  regulons_informed[[regulator]] <- binder(proxy_regulon, m_abscessus_coexpression, m_abscessus_operons, mu_zeta=-3, mu_tau=c(3,3), is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
}

final_regulons_informed <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_informed[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_informed[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_informed[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0){
    print(n)
    final_regulons_informed <- rbind( final_regulons_informed,
    data.frame(
     regulons_informed[[n]]$regulator,
     regulons_informed[[n]]$target_candidate,
     as.numeric(regulons_informed[[n]]$ME),
     as.numeric(regulons_informed[[n]]$PE),
     as.numeric(regulons_informed[[n]]$CM),
     as.numeric(regulons_informed[[n]]$CP),
     as.numeric(regulons_informed[[n]]$mean_theta),
     as.numeric(regulons_informed[[n]]$theta_interval[,1]),
     as.numeric(regulons_informed[[n]]$theta_interval[,2]),
     as.numeric(regulons_informed[[n]]$theta_interval[,3]),
     as.numeric(regulons_informed[[n]]$mean_gamma),
     as.numeric(regulons_informed[[n]]$gamma_interval[,1]),
     as.numeric(regulons_informed[[n]]$gamma_interval[,2]),
     as.numeric(regulons_informed[[n]]$mean_raw_gamma),
     as.numeric(regulons_informed[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_informed[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_informed[[n]]$mean_phi[1]),
     as.numeric(regulons_informed[[n]]$phi_interval[1]),
     as.numeric(regulons_informed[[n]]$phi_interval[2]),
     as.numeric(regulons_informed[[n]]$mean_zeta[1]),
     as.numeric(regulons_informed[[n]]$zeta_interval[1]),
     as.numeric(regulons_informed[[n]]$zeta_interval[2]),
     as.numeric(regulons_informed[[n]]$mean_tau[1]),
     as.numeric(regulons_informed[[n]]$tau_interval[1,1]),
     as.numeric(regulons_informed[[n]]$tau_interval[1,2]),
     as.numeric(regulons_informed[[n]]$mean_tau[2]),
     as.numeric(regulons_informed[[n]]$tau_interval[2,1]),
     as.numeric(regulons_informed[[n]]$tau_interval[2,2]),
     as.numeric(regulons_informed[[n]]$mean_psi[1]),
     as.numeric(regulons_informed[[n]]$psi_interval[1,1]),
     as.numeric(regulons_informed[[n]]$psi_interval[1,2]),
     as.numeric(regulons_informed[[n]]$mean_psi[2]),
     as.numeric(regulons_informed[[n]]$psi_interval[2,1]),
     as.numeric(regulons_informed[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_informed) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_informed, file="RESULTS/regulons_results/final_regulons.informed.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

##################################################
##################################################

###
# Run BINDER.
###

n_chains <- 3
n_iter <- 150
n_warmup <- 100
seed <- 1

#specific:
regulons_specific <- list()
for(n in 1:n_regulators){
  print(n)
  regulator <- all_regulators[n]
  print(regulator)
  proxy_regulon <- all_proxy_regulons[all_proxy_regulons$regulator == regulator, ]
  regulons_specific[[regulator]] <- binder(proxy_regulon, m_abscessus_coexpression, m_abscessus_operons, mu_zeta=-3, sigma_zeta=0.1, mu_tau=c(3,3), sigma_tau=c(0.1,0.1), mu_phi=0, sigma_phi=1, mu_psi=c(0,0), sigma_psi=c(1,1), is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
}


final_regulons_specific <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_specific[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_specific[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_specific[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0){
    print(n)
    final_regulons_specific <- rbind( final_regulons_specific,
    data.frame(
     regulons_specific[[n]]$regulator,
     regulons_specific[[n]]$target_candidate,
     as.numeric(regulons_specific[[n]]$ME),
     as.numeric(regulons_specific[[n]]$PE),
     as.numeric(regulons_specific[[n]]$CM),
     as.numeric(regulons_specific[[n]]$CP),
     as.numeric(regulons_specific[[n]]$mean_theta),
     as.numeric(regulons_specific[[n]]$theta_interval[,1]),
     as.numeric(regulons_specific[[n]]$theta_interval[,2]),
     as.numeric(regulons_specific[[n]]$theta_interval[,3]),
     as.numeric(regulons_specific[[n]]$mean_gamma),
     as.numeric(regulons_specific[[n]]$gamma_interval[,1]),
     as.numeric(regulons_specific[[n]]$gamma_interval[,2]),
     as.numeric(regulons_specific[[n]]$mean_raw_gamma),
     as.numeric(regulons_specific[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_specific[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_specific[[n]]$mean_phi[1]),
     as.numeric(regulons_specific[[n]]$phi_interval[1]),
     as.numeric(regulons_specific[[n]]$phi_interval[2]),
     as.numeric(regulons_specific[[n]]$mean_zeta[1]),
     as.numeric(regulons_specific[[n]]$zeta_interval[1]),
     as.numeric(regulons_specific[[n]]$zeta_interval[2]),
     as.numeric(regulons_specific[[n]]$mean_tau[1]),
     as.numeric(regulons_specific[[n]]$tau_interval[1,1]),
     as.numeric(regulons_specific[[n]]$tau_interval[1,2]),
     as.numeric(regulons_specific[[n]]$mean_tau[2]),
     as.numeric(regulons_specific[[n]]$tau_interval[2,1]),
     as.numeric(regulons_specific[[n]]$tau_interval[2,2]),
     as.numeric(regulons_specific[[n]]$mean_psi[1]),
     as.numeric(regulons_specific[[n]]$psi_interval[1,1]),
     as.numeric(regulons_specific[[n]]$psi_interval[1,2]),
     as.numeric(regulons_specific[[n]]$mean_psi[2]),
     as.numeric(regulons_specific[[n]]$psi_interval[2,1]),
     as.numeric(regulons_specific[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_specific) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_specific, file="RESULTS/regulons_results/final_regulons.specific.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

##################################################
##################################################

final_regulons_vague_common <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_vague[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_vague[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_vague[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
     rstan::get_num_divergent(regulons_informed[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_informed[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_informed[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
     rstan::get_num_divergent(regulons_specific[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_specific[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_specific[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0){
    print(n)
    final_regulons_vague_common <- rbind( final_regulons_vague_common,
    data.frame(
     regulons_vague[[n]]$regulator,
     regulons_vague[[n]]$target_candidate,
     as.numeric(regulons_vague[[n]]$ME),
     as.numeric(regulons_vague[[n]]$PE),
     as.numeric(regulons_vague[[n]]$CM),
     as.numeric(regulons_vague[[n]]$CP),
     as.numeric(regulons_vague[[n]]$mean_theta),
     as.numeric(regulons_vague[[n]]$theta_interval[,1]),
     as.numeric(regulons_vague[[n]]$theta_interval[,2]),
     as.numeric(regulons_vague[[n]]$theta_interval[,3]),
     as.numeric(regulons_vague[[n]]$mean_gamma),
     as.numeric(regulons_vague[[n]]$gamma_interval[,1]),
     as.numeric(regulons_vague[[n]]$gamma_interval[,2]),
     as.numeric(regulons_vague[[n]]$mean_raw_gamma),
     as.numeric(regulons_vague[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_vague[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_vague[[n]]$mean_phi[1]),
     as.numeric(regulons_vague[[n]]$phi_interval[1]),
     as.numeric(regulons_vague[[n]]$phi_interval[2]),
     as.numeric(regulons_vague[[n]]$mean_zeta[1]),
     as.numeric(regulons_vague[[n]]$zeta_interval[1]),
     as.numeric(regulons_vague[[n]]$zeta_interval[2]),
     as.numeric(regulons_vague[[n]]$mean_tau[1]),
     as.numeric(regulons_vague[[n]]$tau_interval[1,1]),
     as.numeric(regulons_vague[[n]]$tau_interval[1,2]),
     as.numeric(regulons_vague[[n]]$mean_tau[2]),
     as.numeric(regulons_vague[[n]]$tau_interval[2,1]),
     as.numeric(regulons_vague[[n]]$tau_interval[2,2]),
     as.numeric(regulons_vague[[n]]$mean_psi[1]),
     as.numeric(regulons_vague[[n]]$psi_interval[1,1]),
     as.numeric(regulons_vague[[n]]$psi_interval[1,2]),
     as.numeric(regulons_vague[[n]]$mean_psi[2]),
     as.numeric(regulons_vague[[n]]$psi_interval[2,1]),
     as.numeric(regulons_vague[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_vague_common) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_vague_common, file="RESULTS/regulons_results/final_regulons.vague.common.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


final_regulons_informed_common <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_vague[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_vague[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_vague[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
    rstan::get_num_divergent(regulons_informed[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_informed[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_informed[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
    rstan::get_num_divergent(regulons_specific[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_specific[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_specific[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0
  ){
    print(n)
    final_regulons_informed_common <- rbind( final_regulons_informed_common,
    data.frame(
     regulons_informed[[n]]$regulator,
     regulons_informed[[n]]$target_candidate,
     as.numeric(regulons_informed[[n]]$ME),
     as.numeric(regulons_informed[[n]]$PE),
     as.numeric(regulons_informed[[n]]$CM),
     as.numeric(regulons_informed[[n]]$CP),
     as.numeric(regulons_informed[[n]]$mean_theta),
     as.numeric(regulons_informed[[n]]$theta_interval[,1]),
     as.numeric(regulons_informed[[n]]$theta_interval[,2]),
     as.numeric(regulons_informed[[n]]$theta_interval[,3]),
     as.numeric(regulons_informed[[n]]$mean_gamma),
     as.numeric(regulons_informed[[n]]$gamma_interval[,1]),
     as.numeric(regulons_informed[[n]]$gamma_interval[,2]),
     as.numeric(regulons_informed[[n]]$mean_raw_gamma),
     as.numeric(regulons_informed[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_informed[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_informed[[n]]$mean_phi[1]),
     as.numeric(regulons_informed[[n]]$phi_interval[1]),
     as.numeric(regulons_informed[[n]]$phi_interval[2]),
     as.numeric(regulons_informed[[n]]$mean_zeta[1]),
     as.numeric(regulons_informed[[n]]$zeta_interval[1]),
     as.numeric(regulons_informed[[n]]$zeta_interval[2]),
     as.numeric(regulons_informed[[n]]$mean_tau[1]),
     as.numeric(regulons_informed[[n]]$tau_interval[1,1]),
     as.numeric(regulons_informed[[n]]$tau_interval[1,2]),
     as.numeric(regulons_informed[[n]]$mean_tau[2]),
     as.numeric(regulons_informed[[n]]$tau_interval[2,1]),
     as.numeric(regulons_informed[[n]]$tau_interval[2,2]),
     as.numeric(regulons_informed[[n]]$mean_psi[1]),
     as.numeric(regulons_informed[[n]]$psi_interval[1,1]),
     as.numeric(regulons_informed[[n]]$psi_interval[1,2]),
     as.numeric(regulons_informed[[n]]$mean_psi[2]),
     as.numeric(regulons_informed[[n]]$psi_interval[2,1]),
     as.numeric(regulons_informed[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_informed_common) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_informed_common, file="RESULTS/regulons_results/final_regulons.informed.common.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


final_regulons_specific_common <- data.frame( stringsAsFactors=FALSE )
for(n in 1:n_regulators){
  if(rstan::get_num_divergent(regulons_vague[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_vague[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_vague[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
    rstan::get_num_divergent(regulons_informed[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_informed[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_informed[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0 &
    rstan::get_num_divergent(regulons_specific[[n]]$model_object) == 0 & length(rstan::get_low_bfmi_chains(regulons_specific[[n]]$model_object)) == 0 & length(which(rstan::summary(regulons_specific[[n]]$model_object)$rstan::summary[, "Rhat"] > 2)) == 0
  ){
    print(n)
    final_regulons_specific_common <- rbind( final_regulons_specific_common,
    data.frame(
    regulons_specific[[n]]$regulator,
    regulons_specific[[n]]$target_candidate,
     as.numeric(regulons_specific[[n]]$ME),
     as.numeric(regulons_specific[[n]]$PE),
     as.numeric(regulons_specific[[n]]$CM),
     as.numeric(regulons_specific[[n]]$CP),
     as.numeric(regulons_specific[[n]]$mean_theta),
     as.numeric(regulons_specific[[n]]$theta_interval[,1]),
     as.numeric(regulons_specific[[n]]$theta_interval[,2]),
     as.numeric(regulons_specific[[n]]$theta_interval[,3]),
     as.numeric(regulons_specific[[n]]$mean_gamma),
     as.numeric(regulons_specific[[n]]$gamma_interval[,1]),
     as.numeric(regulons_specific[[n]]$gamma_interval[,2]),
     as.numeric(regulons_specific[[n]]$mean_raw_gamma),
     as.numeric(regulons_specific[[n]]$raw_gamma_interval[,1]),
     as.numeric(regulons_specific[[n]]$raw_gamma_interval[,2]),
     as.numeric(regulons_specific[[n]]$mean_phi[1]),
     as.numeric(regulons_specific[[n]]$phi_interval[1]),
     as.numeric(regulons_specific[[n]]$phi_interval[2]),
     as.numeric(regulons_specific[[n]]$mean_zeta[1]),
     as.numeric(regulons_specific[[n]]$zeta_interval[1]),
     as.numeric(regulons_specific[[n]]$zeta_interval[2]),
     as.numeric(regulons_specific[[n]]$mean_tau[1]),
     as.numeric(regulons_specific[[n]]$tau_interval[1,1]),
     as.numeric(regulons_specific[[n]]$tau_interval[1,2]),
     as.numeric(regulons_specific[[n]]$mean_tau[2]),
     as.numeric(regulons_specific[[n]]$tau_interval[2,1]),
     as.numeric(regulons_specific[[n]]$tau_interval[2,2]),
     as.numeric(regulons_specific[[n]]$mean_psi[1]),
     as.numeric(regulons_specific[[n]]$psi_interval[1,1]),
     as.numeric(regulons_specific[[n]]$psi_interval[1,2]),
     as.numeric(regulons_specific[[n]]$mean_psi[2]),
     as.numeric(regulons_specific[[n]]$psi_interval[2,1]),
     as.numeric(regulons_specific[[n]]$psi_interval[2,2])
    )
    )
  }
}
colnames(final_regulons_specific_common) <- c("Regulator", "TargetCandidate", "ME", "PE", "CM", "CP", "MeanTheta", "0.025ThetaInterval", "0.5ThetaInterval", "0.975ThetaInterval", "MeanGamma", "0.025GammaInterval", "0.975GammaInterval", "MeanRawGamma", "0.025RawGammaInterval", "0.975RawGammaInterval", "MeanPhi", "0.025PhiInterval", "0.975PhiInterval", "MeanZeta", "0.025ZetaInterval", "0.975ZetaInterval", "MeanTau[1]", "0.025TauInterval[1]", "0.975TauInterval[1]", "MeanTau[2]", "0.025TauInterval[2]", "0.975TauInterval[2]", "MeanPsi[1]", "0.025PsiInterval[1]", "0.975PsiInterval[1]", "MeanPsi[2]", "0.025PsiInterval[2]", "0.975PsiInterval[2]")

write.table(final_regulons_specific_common, file="RESULTS/regulons_results/final_regulons.specific.common.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
