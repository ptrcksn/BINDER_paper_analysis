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
library(rstan)
library(iRafNet)
library(doParallel)
library(ggplot2)
library(plotROC)
library(grid)
library(gridExtra)

#########################
#########################

###
# Parallelize Stan.
###

no_cores <- min(3, (parallel::detectCores()-1))
rstan_options(auto_write=TRUE)
options(mc.cores=(no_cores))

#########################
#########################

###
# DATA.
###

# E. coli:

e_coli_expression <- read.delim("DATA/e_coli/expression.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)

e_coli_locus_tags <- e_coli_expression[, 1]
e_coli_gene_names <- e_coli_expression[, 2]

e_coli_expression <- e_coli_expression[, 4:ncol(e_coli_expression)] # Only keep expression columns.
rownames(e_coli_expression) <- e_coli_locus_tags

no_na_e_coli <- c()
for(m in 1:ncol(e_coli_expression)){if(length(which(is.na(e_coli_expression[, m]))) == 0){no_na_e_coli <- c(no_na_e_coli, m)}}

e_coli_expression_all <- e_coli_expression
e_coli_expression <- e_coli_expression[, no_na_e_coli]



e_coli_fur_proxy_regulon <- read.delim("DATA/e_coli/p_aeruginosa/ortholog_ME_PE_fur.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_fur_regulator <- unique(e_coli_fur_proxy_regulon$regulator)

e_coli_lexA_proxy_regulon <- read.delim("DATA/e_coli/p_aeruginosa/ortholog_ME_PE_lexA.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_lexA_regulator <- unique(e_coli_lexA_proxy_regulon$regulator)

keep_e_coli <- sort(intersect(rownames(e_coli_expression), intersect(e_coli_fur_proxy_regulon$target_candidate, e_coli_lexA_proxy_regulon$target_candidate)))

e_coli_expression_all <- e_coli_expression_all[keep_e_coli, ]
e_coli_expression <- e_coli_expression[keep_e_coli, ]
e_coli_standardised_expression <- apply(t(e_coli_expression), 2, function(x){(x-mean(x))/sd(x)})

e_coli_fur_proxy_regulon <- e_coli_fur_proxy_regulon[e_coli_fur_proxy_regulon$target_candidate %in% keep_e_coli, ]
e_coli_fur_proxy_regulon <- e_coli_fur_proxy_regulon[order(keep_e_coli), ]

e_coli_lexA_proxy_regulon <- e_coli_lexA_proxy_regulon[e_coli_lexA_proxy_regulon$target_candidate %in% keep_e_coli, ]
e_coli_lexA_proxy_regulon <- e_coli_lexA_proxy_regulon[order(keep_e_coli), ]



validated_e_coli_interactions <- read.table("DATA/e_coli/network_tf_gene.txt", header=FALSE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_fur_interactions <- unique(e_coli_fur_proxy_regulon$target_candidate[e_coli_fur_proxy_regulon$target_candidate %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions$V2[validated_e_coli_interactions$V1 == "Fur" & validated_e_coli_interactions$V5 == "Strong"], function(x){which(e_coli_gene_names == x)}))]])
e_coli_lexA_interactions <- unique(e_coli_lexA_proxy_regulon$target_candidate[e_coli_lexA_proxy_regulon$target_candidate %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions$V2[validated_e_coli_interactions$V1 == "LexA" & validated_e_coli_interactions$V5 == "Strong"], function(x){which(e_coli_gene_names == x)}))]])



W <- matrix(rep(0, length(keep_e_coli)^2), nrow=length(keep_e_coli))
rownames(W) <- colnames(W) <- keep_e_coli

e_coli_fur_idx <- which(colnames(e_coli_standardised_expression) == e_coli_fur_regulator)
W[e_coli_fur_regulator, e_coli_fur_proxy_regulon$target_candidate] <- ifelse(e_coli_fur_proxy_regulon$PE == 1 | e_coli_fur_proxy_regulon$ME == 1, 1, 0)
W[e_coli_fur_regulator, e_coli_fur_regulator] <- -sum(W[e_coli_fur_regulator, -e_coli_fur_idx])

e_coli_lexA_idx <- which(colnames(e_coli_standardised_expression) == e_coli_lexA_regulator)
W[e_coli_lexA_regulator, e_coli_lexA_proxy_regulon$target_candidate] <- ifelse(e_coli_lexA_proxy_regulon$PE == 1 | e_coli_lexA_proxy_regulon$ME == 1, 1, 0)
W[e_coli_lexA_regulator, e_coli_lexA_regulator] <- -sum(W[e_coli_lexA_regulator, -e_coli_lexA_idx])

W <- exp(W)



###
# Run iRafNet.
###

cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
e_coli_results <- foreach(i=1:length(keep_e_coli), .combine=rbind) %dopar% if(i != e_coli_fur_idx & i != e_coli_lexA_idx) iRafNet(e_coli_standardised_expression[, c(e_coli_fur_idx, e_coli_lexA_idx, i)], W[c(e_coli_fur_idx, e_coli_lexA_idx, i), c(e_coli_fur_idx, e_coli_lexA_idx, i)], mtry=round(sqrt(length(keep_e_coli)-1)), ntree=1000, colnames(e_coli_standardised_expression)[c(e_coli_fur_idx, e_coli_lexA_idx, i)])
stopCluster(cl)



e_coli_fur_scores <- e_coli_results[e_coli_results$gene1 == e_coli_fur_regulator & e_coli_results$gene2 != e_coli_fur_regulator & e_coli_results$gene2 != e_coli_lexA_regulator, ]
e_coli_fur_frame <- data.frame(Organism="E. coli", Regulator="Fur", RegulatorLocus=e_coli_fur_scores$gene1, TargetCandidate=e_coli_fur_scores$gene2, Score=e_coli_fur_scores$importance, Interaction=ifelse(e_coli_fur_scores$gene2 %in% e_coli_fur_interactions, 1, 0), Model="iRafNet", stringsAsFactors=FALSE)

e_coli_lexA_scores <- e_coli_results[e_coli_results$gene1 == e_coli_lexA_regulator & e_coli_results$gene2 != e_coli_fur_regulator & e_coli_results$gene2 != e_coli_lexA_regulator, ]
e_coli_lexA_frame <- data.frame(Organism="E. coli", Regulator="LexA", RegulatorLocus=e_coli_lexA_scores$gene1, TargetCandidate=e_coli_lexA_scores$gene2, Score=e_coli_lexA_scores$importance, Interaction=ifelse(e_coli_lexA_scores$gene2 %in% e_coli_lexA_interactions, 1, 0), Model="iRafNet", stringsAsFactors=FALSE)



###
# DATA.
###

# B. subtilis:

b_subtilis_expression <- read.delim("DATA/b_subtilis/expression.csv", sep=",", header=TRUE, row.names=2, fill=TRUE, stringsAsFactors=FALSE)

b_subtilis_expression <- b_subtilis_expression[, 13:dim(b_subtilis_expression)[2]] # Only keep expression columns.

no_na_b_subtilis <- c()
for(m in 1:ncol(b_subtilis_expression)){if(length(which(is.na(b_subtilis_expression[, m]))) == 0){no_na_b_subtilis <- c(no_na_b_subtilis, m)}}

b_subtilis_expression <- b_subtilis_expression[, no_na_b_subtilis]



b_subtilis_fur_proxy_regulon <- read.delim("DATA/b_subtilis/l_monocytogenes/ortholog_ME_PE_fur.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_fur_regulator <- unique(b_subtilis_fur_proxy_regulon$regulator)

b_subtilis_lexA_proxy_regulon <- read.delim("DATA/b_subtilis/l_monocytogenes/ortholog_ME_PE_lexA.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_lexA_regulator <- unique(b_subtilis_lexA_proxy_regulon$regulator)

keep_b_subtilis <- sort(intersect(rownames(b_subtilis_expression), intersect(b_subtilis_fur_proxy_regulon$target_candidate, b_subtilis_lexA_proxy_regulon$target_candidate)))

b_subtilis_expression <- b_subtilis_expression[keep_b_subtilis, ]
b_subtilis_standardised_expression <- apply(t(b_subtilis_expression), 2, function(x){(x-mean(x))/sd(x)})

b_subtilis_fur_proxy_regulon <- b_subtilis_fur_proxy_regulon[b_subtilis_fur_proxy_regulon$target_candidate %in% keep_b_subtilis, ]
b_subtilis_fur_proxy_regulon <- b_subtilis_fur_proxy_regulon[order(keep_b_subtilis), ]

b_subtilis_lexA_proxy_regulon <- b_subtilis_lexA_proxy_regulon[b_subtilis_lexA_proxy_regulon$target_candidate %in% keep_b_subtilis, ]
b_subtilis_lexA_proxy_regulon <- b_subtilis_lexA_proxy_regulon[order(keep_b_subtilis), ]



validated_b_subtilis_interactions <- read.table("DATA/b_subtilis/regulations.csv", header=TRUE, sep=",", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_fur_interactions <- unique(validated_b_subtilis_interactions$locus.tag[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_fur_regulator])
b_subtilis_lexA_interactions <- unique(validated_b_subtilis_interactions$locus.tag[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_lexA_regulator])



W <- matrix(rep(0, length(keep_b_subtilis)^2), nrow=length(keep_b_subtilis))
rownames(W) <- colnames(W) <- keep_b_subtilis

b_subtilis_fur_idx <- which(colnames(b_subtilis_standardised_expression) == b_subtilis_fur_regulator)
W[b_subtilis_fur_regulator, b_subtilis_fur_proxy_regulon$target_candidate] <- ifelse(b_subtilis_fur_proxy_regulon$PE == 1 | b_subtilis_fur_proxy_regulon$ME == 1, 1, 0)
W[b_subtilis_fur_regulator, b_subtilis_fur_regulator] <- -sum(W[b_subtilis_fur_regulator, -b_subtilis_fur_idx])

b_subtilis_lexA_idx <- which(colnames(b_subtilis_standardised_expression) == b_subtilis_lexA_regulator)
W[b_subtilis_lexA_regulator, b_subtilis_lexA_proxy_regulon$target_candidate] <- ifelse(b_subtilis_lexA_proxy_regulon$PE == 1 | b_subtilis_lexA_proxy_regulon$ME == 1, 1, 0)
W[b_subtilis_lexA_regulator, b_subtilis_lexA_regulator] <- -sum(W[b_subtilis_lexA_regulator, -b_subtilis_lexA_idx])

W <- exp(W)



###
# Run iRafNet.
###

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
b_subtilis_results <- foreach(i=1:length(keep_b_subtilis), .combine=rbind) %dopar% if(i != b_subtilis_fur_idx & i != b_subtilis_lexA_idx) iRafNet(b_subtilis_standardised_expression[, c(b_subtilis_fur_idx, b_subtilis_lexA_idx, i)], W[c(b_subtilis_fur_idx, b_subtilis_lexA_idx, i), c(b_subtilis_fur_idx, b_subtilis_lexA_idx, i)], mtry=round(sqrt(length(keep_b_subtilis)-1)), ntree=1000, colnames(b_subtilis_standardised_expression)[c(b_subtilis_fur_idx, b_subtilis_lexA_idx, i)])
stopCluster(cl)



b_subtilis_fur_scores <- b_subtilis_results[b_subtilis_results$gene1 == b_subtilis_fur_regulator & b_subtilis_results$gene2 != b_subtilis_fur_regulator & b_subtilis_results$gene2 != b_subtilis_lexA_regulator, ]
b_subtilis_fur_frame <- data.frame(Organism="B. subtilis", Regulator="Fur", RegulatorLocus=b_subtilis_fur_scores$gene1, TargetCandidate=b_subtilis_fur_scores$gene2, Score=b_subtilis_fur_scores$importance, Interaction=ifelse(b_subtilis_fur_scores$gene2 %in% b_subtilis_fur_interactions, 1, 0), Model="iRafNet", stringsAsFactors=FALSE)

b_subtilis_lexA_scores <- b_subtilis_results[b_subtilis_results$gene1 == b_subtilis_lexA_regulator & b_subtilis_results$gene2 != b_subtilis_fur_regulator & b_subtilis_results$gene2 != b_subtilis_lexA_regulator, ]
b_subtilis_lexA_frame <- data.frame(Organism="B. subtilis", Regulator="LexA", RegulatorLocus=b_subtilis_lexA_scores$gene1, TargetCandidate=b_subtilis_lexA_scores$gene2, Score=b_subtilis_lexA_scores$importance, Interaction=ifelse(b_subtilis_lexA_scores$gene2 %in% b_subtilis_lexA_interactions, 1, 0), Model="iRafNet", stringsAsFactors=FALSE)





###
# DATA.
###

regulon_data <- list()

# E. coli:

e_coli_coexpression_all <- compute_coexpression(e_coli_expression_all)
e_coli_coexpression <- compute_coexpression(e_coli_expression)

regulon_data[["e_coli_fur"]]$organism <- "E. coli"
regulon_data[["e_coli_fur"]]$regulator <- "Fur"
regulon_data[["e_coli_fur"]]$regulator_locus <- e_coli_fur_regulator
regulon_data[["e_coli_fur"]]$proxy_regulon <- e_coli_fur_proxy_regulon
regulon_data[["e_coli_fur"]]$coexpression_all <- e_coli_coexpression_all
regulon_data[["e_coli_fur"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_fur"]]$interactions <- e_coli_fur_interactions

regulon_data[["e_coli_lexA"]]$organism <- "E. coli"
regulon_data[["e_coli_lexA"]]$regulator <- "LexA"
regulon_data[["e_coli_lexA"]]$regulator_locus <- e_coli_lexA_regulator
regulon_data[["e_coli_lexA"]]$proxy_regulon <- e_coli_lexA_proxy_regulon
regulon_data[["e_coli_lexA"]]$coexpression_all <- e_coli_coexpression_all
regulon_data[["e_coli_lexA"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_lexA"]]$interactions <- e_coli_lexA_interactions


# B. subtilis:

b_subtilis_coexpression <- compute_coexpression(b_subtilis_expression)

regulon_data[["b_subtilis_fur"]]$organism <- "B. subtilis"
regulon_data[["b_subtilis_fur"]]$regulator <- "Fur"
regulon_data[["b_subtilis_fur"]]$regulator_locus <- b_subtilis_fur_regulator
regulon_data[["b_subtilis_fur"]]$proxy_regulon <- b_subtilis_fur_proxy_regulon
regulon_data[["b_subtilis_fur"]]$coexpression <- b_subtilis_coexpression
regulon_data[["b_subtilis_fur"]]$interactions <- b_subtilis_fur_interactions

regulon_data[["b_subtilis_lexA"]]$organism <- "B. subtilis"
regulon_data[["b_subtilis_lexA"]]$regulator <- "LexA"
regulon_data[["b_subtilis_lexA"]]$regulator_locus <- b_subtilis_lexA_regulator
regulon_data[["b_subtilis_lexA"]]$expression <- b_subtilis_expression
regulon_data[["b_subtilis_lexA"]]$proxy_regulon <- b_subtilis_lexA_proxy_regulon
regulon_data[["b_subtilis_lexA"]]$coexpression <- b_subtilis_coexpression
regulon_data[["b_subtilis_lexA"]]$interactions <- b_subtilis_lexA_interactions



###
# BINDER.
###

all_df_binder_all <- c()
all_df_binder_muphi10 <- c()
all_df_binder <- c()
all_df_non_auxiliary <- c()
all_df_regression <- c()
n_regulons <- length(regulon_data)
for(n in 1:n_regulons){

  organism <- regulon_data[[n]]$organism
  regulator <- regulon_data[[n]]$regulator
  regulator_locus <- regulon_data[[n]]$regulator_locus
  proxy_regulon <- regulon_data[[n]]$proxy_regulon
  coexpression <- regulon_data[[n]]$coexpression
  interactions <- regulon_data[[n]]$interactions

  print(regulator)

  if(organism == "E. coli"){
    coexpression_all <-  regulon_data[[n]]$coexpression_all

    results_binder_all <- binder(proxy_regulon, coexpression_all, is_coexpression=TRUE,  mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_phi=1, sigma_phi=0.1, mu_psi=c(0,0), sigma_psi=c(3,3), chains=3, seed=1)[[1]]
    df_binder_all <- data.frame(Organism=organism, Regulator=regulator, RegulatorLocus=regulator_locus, TargetCandidate=results_binder_all$target_candidate, Score=results_binder_all$theta_interval[, 2], Interaction=ifelse(results_binder_all$target_candidate %in% interactions, 1, 0), Model="BINDER (all)", ME=results_binder_all$ME, PE=results_binder_all$PE, CM=results_binder_all$CM, CP=results_binder_all$CP, stringsAsFactors=FALSE)
    all_df_binder_all <- rbind(all_df_binder_all, df_binder_all)

    if(regulator == "Fur"){
      results_binder_muphi10 <- binder(proxy_regulon, coexpression, is_coexpression=TRUE,  mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_phi=10, sigma_phi=0.1, mu_psi=c(0,0), sigma_psi=c(3,3), chains=3, seed=1)[[1]]
      df_binder_muphi10 <- data.frame(Organism=organism, Regulator=regulator, RegulatorLocus=regulator_locus, TargetCandidate=results_binder_muphi10$target_candidate, Score=results_binder_muphi10$theta_interval[, 2], Interaction=ifelse(results_binder_muphi10$target_candidate %in% interactions, 1, 0), Model="BINDER [Informative p(phi)]", ME=results_binder_muphi10$ME, PE=results_binder_muphi10$PE, CM=results_binder_muphi10$CM, CP=results_binder_muphi10$CP, stringsAsFactors=FALSE)
      all_df_binder_muphi10 <- rbind(all_df_binder_muphi10, df_binder_muphi10)
    }
  }

  results_binder <- binder(proxy_regulon, coexpression, is_coexpression=TRUE, mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_phi=1, sigma_phi=0.1, mu_psi=c(0,0), sigma_psi=c(3,3), chains=3, seed=1)[[1]]
  df_binder <- data.frame(Organism=organism, Regulator=regulator, RegulatorLocus=regulator_locus, TargetCandidate=results_binder$target_candidate, Score=results_binder$theta_interval[, 2], Interaction=ifelse(results_binder$target_candidate %in% interactions, 1, 0), Model="BINDER", ME=results_binder$ME, PE=results_binder$PE, CM=results_binder$CM, CP=results_binder$CP, stringsAsFactors=FALSE)
  all_df_binder <- rbind(all_df_binder, df_binder)

  results_non_auxiliary <- binder(proxy_regulon, coexpression, is_coexpression=TRUE, mu_psi=c(0,0), sigma_psi=c(3,3), model="non_auxiliary", chains=3, seed=1)[[1]]
  df_non_auxiliary <- data.frame(Organism=organism, Regulator=regulator, RegulatorLocus=regulator_locus, TargetCandidate=results_non_auxiliary$target_candidate, Score=results_non_auxiliary$theta_interval[, 2], Interaction=ifelse(results_non_auxiliary$target_candidate %in% interactions, 1, 0), Model="Non-auxiliary", ME=results_non_auxiliary$ME, PE=results_non_auxiliary$PE, CM=results_non_auxiliary$CM, CP=results_non_auxiliary$CP, stringsAsFactors=FALSE)
  all_df_non_auxiliary <- rbind(all_df_non_auxiliary, df_non_auxiliary)

  results_regression <- binder(proxy_regulon, coexpression, is_coexpression=TRUE, mu_zeta=0, sigma_zeta=3, mu_tau=c(0,0), sigma_tau=c(3,3), mu_psi=c(0,0), sigma_psi=c(3,3), model="deterministic", chains=3, seed=1)[[1]]
  df_regression <- data.frame(Organism=organism, Regulator=regulator, RegulatorLocus=regulator_locus, TargetCandidate=results_regression$target_candidate, Score=results_regression$theta_interval[, 2], Interaction=ifelse(results_regression$target_candidate %in% interactions, 1, 0), Model="Deterministic", ME=results_regression$ME, PE=results_regression$PE, CM=results_regression$CM, CP=results_regression$CP, stringsAsFactors=FALSE)
  all_df_regression <- rbind(all_df_regression, df_regression)

}



all_df <- rbind(all_df_binder_all[, 1:7], all_df_binder_muphi10[1:7], all_df_binder[, 1:7], all_df_non_auxiliary[, 1:7], all_df_regression[, 1:7], e_coli_fur_frame, e_coli_lexA_frame, b_subtilis_fur_frame, b_subtilis_lexA_frame)
all_df$combination <- paste(all_df$Regulator, all_df$Organism, all_df$Model, sep="-")

base_palette <- c("#54b8d3", "#EFC69B", "#CE8964", "#832232", "grey", "#2E0219")

combinations <- unique(all_df$combination)
n_combinations <- length(combinations)
for(n in 1:n_combinations){
  current_combination <- combinations[n]
  print(current_combination)
  print("AUC")
  print(calc_auc(ggplot(all_df[all_df$combination == current_combination, ], aes(d=Interaction, m=Score, color=combination)) + geom_roc(n.cuts=0))[,"AUC"])
  print("Interactions Mean Score")
  print(mean(all_df[all_df$combination == current_combination & all_df$Interaction == 1, ]$Score))
  print("Non-interactions Mean Score")
  print(mean(all_df[all_df$combination == current_combination & all_df$Interaction == 0, ]$Score))
  print("----------")
}

ggplot(all_df, aes(d=Interaction, m=Score, color=factor(Model, levels=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)")), size=factor(Model, levels=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)")))) +
  geom_abline(aes(intercept=0, slope=1), linetype=2) +
  geom_roc(n.cuts=0) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  facet_grid(Regulator ~ Organism) +
  theme_minimal() +
  theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=12), legend.text=element_text(size=10), strip.text=element_text(size=12), legend.position="bottom") +
  scale_color_manual(name="Model", breaks=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)"), labels=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", expression(paste("BINDER [Informative ",p(phi),"]")), "BINDER (all)"), values=base_palette) +
  scale_size_manual(name="Model", breaks=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)"), labels=c("iRafNet", "Deterministic", "Non-auxiliary", "BINDER", expression(paste("BINDER [Informative ",p(phi),"]")), "BINDER (all)"), values=c(0.5, 0.4, 0.4, 0.5, 0.5, 0.5)) +
  coord_fixed(ratio=1) +
  guides(shape=guide_legend(override.aes=list(size=1)))

all_df <- rbind(all_df_binder_all, all_df_binder_muphi10, all_df_binder, all_df_non_auxiliary, all_df_regression)
all_df$Interaction <- ifelse(all_df$Interaction == 1, "Interaction", "No Interaction")
all_df$AuxiliaryEvidence <- ifelse(all_df$PE == 1 | all_df$ME == 1, "Auxiliary Evidence", "No Auxiliary Evidence")

base_palette <- c("#EFC69B", "#CE8964", "#832232", "grey", "#2E0219")

ggplot(all_df) +
  geom_boxplot(aes(x=Interaction, y=Score, color=factor(Model, levels=c("Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)")), fill=factor(Model, levels=c("Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)"))), alpha=0.7) +
  xlab("Interaction") +
  ylab(expression(paste(theta^"50%"))) +
  facet_grid(Regulator ~ Organism) +
  theme_minimal() +
  theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=12), legend.text=element_text(size=10), strip.text=element_text(size=12), axis.text.y=element_text(angle=90), legend.position="bottom") +
  scale_fill_manual(name="Model", breaks=c("Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)"), labels=c("Deterministic", "Non-auxiliary", "BINDER", expression(paste("BINDER [Informative ",p(phi),"]")), "BINDER (all)"), values=base_palette) +
  scale_color_manual(name="Model", breaks=c("Deterministic", "Non-auxiliary", "BINDER", "BINDER [Informative p(phi)]", "BINDER (all)"), labels=c("Deterministic", "Non-auxiliary", "BINDER", expression(paste("BINDER [Informative ",p(phi),"]")), "BINDER (all)"), values=base_palette) +
  guides(shape=guide_legend(override.aes=list(size=1)), fill=guide_legend(nrow=2,byrow=TRUE)) +
  coord_flip()

set.seed(1)

ggplot(all_df[all_df$Organism == "B. subtilis" & all_df$Regulator == "LexA" & all_df$AuxiliaryEvidence == "No Auxiliary Evidence", ]) +
  geom_jitter(aes(x=CM, y=CP, fill=Score), shape=22, size=3, alpha=0.5, width=0.1, height=0.1) +
  xlab("CM") +
  ylab("CP") +
  facet_grid(factor(Model, levels=c("Deterministic", "Non-auxiliary", "BINDER")) ~ Interaction) +
  theme_minimal() +
  theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=12), legend.text=element_text(size=10), strip.text=element_text(size=12)) +
  scale_fill_gradientn(name=expression(paste(theta["lexA"]^"50%")), colours=c("#54b8d3", "#EFC69B", "#CE8964", "#832232", "#2E0219"), limits=c(0,1)) +
  coord_fixed(ratio=1)
