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
library(GENIE3)
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
options(mc.cores=(no_cores))

#########################
#########################

###
# DATA.
###

regulon_data <- list()

# E. coli (all):
e_coli_expression <- read.delim("DATA/e_coli/expression.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)
e_coli_locus_tags <- e_coli_expression[, 1]
e_coli_gene_names <- e_coli_expression[, 2]
e_coli_expression <- e_coli_expression[, 4:ncol(e_coli_expression)] # Only keep expression columns.
rownames(e_coli_expression) <- e_coli_locus_tags

e_coli_operons <- read.delim("DATA/e_coli/operons.opr", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_operons <- e_coli_operons[order(e_coli_operons$Synonym), ]

keep_e_coli <- intersect(e_coli_operons$Synonym, rownames(e_coli_expression))
e_coli_expression <- e_coli_expression[keep_e_coli, ]
e_coli_operons <- e_coli_operons[e_coli_operons$Synonym %in% keep_e_coli, ]
e_coli_coexpression <- compute_coexpression(e_coli_expression, keep_e_coli)

validated_e_coli_interactions <- read.table("DATA/e_coli/network_tf_gene.txt", header=FALSE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

e_coli_fur_proxy_regulon <- read.delim("DATA/e_coli/paeruginosa/ortholog_ME_PE_fur.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_fur_regulator <- unique(e_coli_fur_proxy_regulon$regulator)

validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_fur_regulator)]), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V5 == "Strong", ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[!duplicated(validated_e_coli_interactions_regulator$V2), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V2 %in% unique(e_coli_gene_names[which(e_coli_locus_tags %in% keep_e_coli)]), ]
e_coli_fur_interaction_status <- ifelse(keep_e_coli %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(e_coli_gene_names == x)}))], 1, 0)

regulon_data[["e_coli_fur_all"]]$organism <- "E. coli"
regulon_data[["e_coli_fur_all"]]$regulator_name <- "fur"
regulon_data[["e_coli_fur_all"]]$regulator <- e_coli_fur_regulator
regulon_data[["e_coli_fur_all"]]$expression <- e_coli_expression
regulon_data[["e_coli_fur_all"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_fur_all"]]$operons <- e_coli_operons
regulon_data[["e_coli_fur_all"]]$proxy_regulon <- e_coli_fur_proxy_regulon
regulon_data[["e_coli_fur_all"]]$interaction_status <- e_coli_fur_interaction_status

e_coli_lexA_proxy_regulon <- read.delim("DATA/e_coli/paeruginosa/ortholog_ME_PE_lexA.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_lexA_regulator <- unique(e_coli_lexA_proxy_regulon$regulator)

validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_lexA_regulator)]), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V5 == "Strong", ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[!duplicated(validated_e_coli_interactions_regulator$V2), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V2 %in% unique(e_coli_gene_names[which(e_coli_locus_tags %in% keep_e_coli)]), ]
e_coli_lexA_interaction_status <- ifelse(keep_e_coli %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(e_coli_gene_names == x)}))], 1, 0)

regulon_data[["e_coli_lexA_all"]]$organism <- "E. coli"
regulon_data[["e_coli_lexA_all"]]$regulator_name <- "lexA"
regulon_data[["e_coli_lexA_all"]]$regulator <- e_coli_lexA_regulator
regulon_data[["e_coli_lexA_all"]]$expression <- e_coli_expression
regulon_data[["e_coli_lexA_all"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_lexA_all"]]$operons <- e_coli_operons
regulon_data[["e_coli_lexA_all"]]$proxy_regulon <- e_coli_lexA_proxy_regulon
regulon_data[["e_coli_lexA_all"]]$interaction_status <- e_coli_lexA_interaction_status



# E. coli:
e_coli_expression <- read.delim("DATA/e_coli/expression.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)
e_coli_locus_tags <- e_coli_expression[, 1]
e_coli_gene_names <- e_coli_expression[, 2]
e_coli_expression <- e_coli_expression[, 4:ncol(e_coli_expression)] # Only keep expression columns.
rownames(e_coli_expression) <- e_coli_locus_tags
no_na_e_coli <- c()
for(m in 1:ncol(e_coli_expression)){if(length(which(is.na(e_coli_expression[, m]))) == 0){no_na_e_coli <- c(no_na_e_coli, m)}}
e_coli_expression <- e_coli_expression[, no_na_e_coli]

e_coli_operons <- read.delim("DATA/e_coli/operons.opr", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_operons <- e_coli_operons[order(e_coli_operons$Synonym), ]

keep_e_coli <- intersect(e_coli_operons$Synonym, rownames(e_coli_expression))
e_coli_expression <- e_coli_expression[keep_e_coli, ]
e_coli_operons <- e_coli_operons[e_coli_operons$Synonym %in% keep_e_coli, ]
e_coli_coexpression <- compute_coexpression(e_coli_expression, keep_e_coli)

validated_e_coli_interactions <- read.table("DATA/e_coli/network_tf_gene.txt", header=FALSE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

e_coli_fur_proxy_regulon <- read.delim("DATA/e_coli/paeruginosa/ortholog_ME_PE_fur.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_fur_regulator <- unique(e_coli_fur_proxy_regulon$regulator)

validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_fur_regulator)]), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V5 == "Strong", ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[!duplicated(validated_e_coli_interactions_regulator$V2), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V2 %in% unique(e_coli_gene_names[which(e_coli_locus_tags %in% keep_e_coli)]), ]
e_coli_fur_interaction_status <- ifelse(keep_e_coli %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(e_coli_gene_names == x)}))], 1, 0)

regulon_data[["e_coli_fur"]]$organism <- "E. coli"
regulon_data[["e_coli_fur"]]$regulator_name <- "fur"
regulon_data[["e_coli_fur"]]$regulator <- e_coli_fur_regulator
regulon_data[["e_coli_fur"]]$expression <- e_coli_expression
regulon_data[["e_coli_fur"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_fur"]]$operons <- e_coli_operons
regulon_data[["e_coli_fur"]]$proxy_regulon <- e_coli_fur_proxy_regulon
regulon_data[["e_coli_fur"]]$interaction_status <- e_coli_fur_interaction_status

e_coli_lexA_proxy_regulon <- read.delim("DATA/e_coli/paeruginosa/ortholog_ME_PE_lexA.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
e_coli_lexA_regulator <- unique(e_coli_lexA_proxy_regulon$regulator)

validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_lexA_regulator)]), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V5 == "Strong", ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[!duplicated(validated_e_coli_interactions_regulator$V2), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions_regulator[validated_e_coli_interactions_regulator$V2 %in% unique(e_coli_gene_names[which(e_coli_locus_tags %in% keep_e_coli)]), ]
e_coli_lexA_interaction_status <- ifelse(keep_e_coli %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(e_coli_gene_names == x)}))], 1, 0)

regulon_data[["e_coli_lexA"]]$organism <- "E. coli"
regulon_data[["e_coli_lexA"]]$regulator_name <- "lexA"
regulon_data[["e_coli_lexA"]]$regulator <- e_coli_lexA_regulator
regulon_data[["e_coli_lexA"]]$expression <- e_coli_expression
regulon_data[["e_coli_lexA"]]$coexpression <- e_coli_coexpression
regulon_data[["e_coli_lexA"]]$operons <- e_coli_operons
regulon_data[["e_coli_lexA"]]$proxy_regulon <- e_coli_lexA_proxy_regulon
regulon_data[["e_coli_lexA"]]$interaction_status <- e_coli_lexA_interaction_status



# B. subtilis:
b_subtilis_expression <- read.delim("DATA/b_subtilis/expression.csv", sep=",", header=TRUE, row.names=2, fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_expression <- b_subtilis_expression[, 13:dim(b_subtilis_expression)[2]] # Only keep expression columns.
no_na_b_subtilis <- c()
for(m in 1:ncol(b_subtilis_expression)){if(length(which(is.na(b_subtilis_expression[, m]))) == 0){no_na_b_subtilis <- c(no_na_b_subtilis, m)}}
b_subtilis_expression <- b_subtilis_expression[, no_na_b_subtilis]

b_subtilis_operons <- read.delim("DATA/b_subtilis/operons.opr", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_operons <- b_subtilis_operons[order(b_subtilis_operons$Synonym), ]

keep_b_subtilis <- intersect(b_subtilis_operons$Synonym, rownames(b_subtilis_expression))
b_subtilis_expression <- b_subtilis_expression[keep_b_subtilis, ]
b_subtilis_operons <- b_subtilis_operons[b_subtilis_operons$Synonym %in% keep_b_subtilis, ]
b_subtilis_coexpression <- compute_coexpression(b_subtilis_expression, keep_b_subtilis)

validated_b_subtilis_interactions <- read.table("DATA/b_subtilis/regulations.csv", header=TRUE, sep=",", fill=TRUE, stringsAsFactors=FALSE)

b_subtilis_fur_proxy_regulon <- read.delim("DATA/b_subtilis/l_monocytogenes/ortholog_ME_PE_fur.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_fur_regulator <- unique(b_subtilis_fur_proxy_regulon$regulator)

validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_fur_regulator, ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions_regulator[!duplicated(validated_b_subtilis_interactions_regulator$locus.tag), ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions_regulator[validated_b_subtilis_interactions_regulator$locus.tag %in% keep_b_subtilis, ]
b_subtilis_fur_interaction_status <- ifelse(keep_b_subtilis %in% validated_b_subtilis_interactions_regulator$locus.tag, 1, 0)

regulon_data[["b_subtilis_fur"]]$organism <- "B. subtilis"
regulon_data[["b_subtilis_fur"]]$regulator_name <- "fur"
regulon_data[["b_subtilis_fur"]]$regulator <- b_subtilis_fur_regulator
regulon_data[["b_subtilis_fur"]]$expression <- b_subtilis_expression
regulon_data[["b_subtilis_fur"]]$coexpression <- b_subtilis_coexpression
regulon_data[["b_subtilis_fur"]]$operons <- b_subtilis_operons
regulon_data[["b_subtilis_fur"]]$proxy_regulon <- b_subtilis_fur_proxy_regulon
regulon_data[["b_subtilis_fur"]]$interaction_status <- b_subtilis_fur_interaction_status

b_subtilis_lexA_proxy_regulon <- read.delim("DATA/b_subtilis/l_monocytogenes/ortholog_ME_PE_lexA.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
b_subtilis_lexA_regulator <- unique(b_subtilis_lexA_proxy_regulon$regulator)

validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_lexA_regulator, ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions_regulator[!duplicated(validated_b_subtilis_interactions_regulator$locus.tag), ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions_regulator[validated_b_subtilis_interactions_regulator$locus.tag %in% keep_b_subtilis, ]
b_subtilis_lexA_interaction_status <- ifelse(keep_b_subtilis %in% validated_b_subtilis_interactions_regulator$locus.tag, 1, 0)

regulon_data[["b_subtilis_lexA"]]$organism <- "B. subtilis"
regulon_data[["b_subtilis_lexA"]]$regulator_name <- "lexA"
regulon_data[["b_subtilis_lexA"]]$regulator <- b_subtilis_lexA_regulator
regulon_data[["b_subtilis_lexA"]]$expression <- b_subtilis_expression
regulon_data[["b_subtilis_lexA"]]$coexpression <- b_subtilis_coexpression
regulon_data[["b_subtilis_lexA"]]$operons <- b_subtilis_operons
regulon_data[["b_subtilis_lexA"]]$proxy_regulon <- b_subtilis_lexA_proxy_regulon
regulon_data[["b_subtilis_lexA"]]$interaction_status <- b_subtilis_lexA_interaction_status

#########################
#########################

###
# GENIE3.
###

# E. coli:
genie_e_coli_regulators <- intersect(e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions$V1, function(x){which(tolower(e_coli_gene_names) == tolower(x))}))], keep_e_coli)
genie_e_coli_expression <- as.matrix(e_coli_expression)
weightMat_e_coli <- GENIE3(genie_e_coli_expression, regulators=genie_e_coli_regulators, nTrees=100)
linkList_e_coli <- getLinkList(weightMat_e_coli)
all_df_genie <- c()

# fur:
linkList_b0683 <- linkList_e_coli[linkList_e_coli$regulatoryGene == "b0683", ]
linkList_b0683 <- linkList_b0683[order(linkList_b0683$targetGene), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_fur_regulator)]), ]
interaction_status <- ifelse(linkList_b0683$targetGene %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(tolower(e_coli_gene_names) == tolower(x))}))], 1, 0)
df_genie <- data.frame("Regulator"=linkList_b0683$regulatoryGene, "TargetCandidate"=linkList_b0683$targetGene, "Score"=linkList_b0683$weight, "Interaction"=interaction_status, regulator_name="fur", organism="E. coli", model="GENIE3")
all_df_genie <- rbind(all_df_genie, df_genie)
# lexA:
linkList_b4043 <- linkList_e_coli[linkList_e_coli$regulatoryGene == "b4043", ]
linkList_b4043 <- linkList_b4043[order(linkList_b4043$targetGene), ]
validated_e_coli_interactions_regulator <- validated_e_coli_interactions[tolower(validated_e_coli_interactions$V1) == tolower(e_coli_gene_names[which(e_coli_locus_tags == e_coli_lexA_regulator)]), ]
interaction_status <- ifelse(linkList_b4043$targetGene %in% e_coli_locus_tags[as.numeric(sapply(validated_e_coli_interactions_regulator$V2, function(x){which(tolower(e_coli_gene_names) == tolower(x))}))], 1, 0)
df_genie <- data.frame("Regulator"=linkList_b4043$regulatoryGene, "TargetCandidate"=linkList_b4043$targetGene, "Score"=linkList_b4043$weight, "Interaction"=interaction_status, regulator_name="lexA", organism="E. coli", model="GENIE3")
all_df_genie <- rbind(all_df_genie, df_genie)

# B. subtilis:
genie_b_subtilis_regulators <- intersect(validated_b_subtilis_interactions$regulator.locus.tag, keep_b_subtilis)
genie_b_subtilis_expression <- as.matrix(b_subtilis_expression)
weightMat_b_subtilis <- GENIE3(genie_b_subtilis_expression, regulators=genie_b_subtilis_regulators, nTrees=100)
linkList_b_subtilis <- getLinkList(weightMat_b_subtilis)

# fur:
linkList_BSU23520 <- linkList_b_subtilis[linkList_b_subtilis$regulatoryGene == "BSU23520", ]
linkList_BSU23520 <- linkList_BSU23520[order(linkList_BSU23520$targetGene), ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_fur_regulator, ]
interaction_status <- ifelse(linkList_BSU23520$targetGene %in% validated_b_subtilis_interactions_regulator$locus.tag, 1, 0)
df_genie <- data.frame("Regulator"=linkList_BSU23520$regulatoryGene, "TargetCandidate"=linkList_BSU23520$targetGene, "Score"=linkList_BSU23520$weight, "Interaction"=interaction_status, regulator_name="fur", organism="B. subtilis", model="GENIE3")
all_df_genie <- rbind(all_df_genie, df_genie)
# lexA:
linkList_BSU17850 <- linkList_b_subtilis[linkList_b_subtilis$regulatoryGene == "BSU17850", ]
linkList_BSU17850 <- linkList_BSU17850[order(linkList_BSU17850$targetGene), ]
validated_b_subtilis_interactions_regulator <- validated_b_subtilis_interactions[validated_b_subtilis_interactions$regulator.locus.tag == b_subtilis_lexA_regulator, ]
interaction_status <- ifelse(linkList_BSU17850$targetGene %in% validated_b_subtilis_interactions_regulator$locus.tag, 1, 0)
df_genie <- data.frame("Regulator"=linkList_BSU17850$regulatoryGene, "TargetCandidate"=linkList_BSU17850$targetGene, "Score"=linkList_BSU17850$weight, "Interaction"=interaction_status, regulator_name="lexA", organism="B. subtilis", model="GENIE3")
all_df_genie <- rbind(all_df_genie, df_genie)

#########################
#########################

n_chains <- 3
n_iter <- 150
n_warmup <- 100
seed <- 1

###
# BINDER.
###

all_df_binder <- c()
all_df_non_auxiliary <- c()
all_df_deterministic <- c()
n_regulons <- length(regulon_data)
for(n in 1:n_regulons){
  organism <- regulon_data[[n]]$organism
  regulator_name <- regulon_data[[n]]$regulator_name
  regulator <- regulon_data[[n]]$regulator
  proxy_regulon <- regulon_data[[n]]$proxy_regulon
  expression <- regulon_data[[n]]$expression
  coexpression <- regulon_data[[n]]$coexpression
  operons <- regulon_data[[n]]$operons
  interaction_status <- regulon_data[[n]]$interaction_status

  print(regulator)

  if(n <= 2){
    results_binder <- binder(proxy_regulon, coexpression, operons, is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
    df_binder <- data.frame("Regulator"=regulator, "TargetCandidate"=results_binder$target_candidate, "Score"=results_binder$theta_interval[, 2], "Interaction"=interaction_status, regulator_name=regulator_name, organism=organism, model="BINDER (all)", ME=results_binder$ME, PE=results_binder$PE, CM=results_binder$CM, CP=results_binder$CP)
    all_df_binder <- rbind(all_df_binder, df_binder)
  }else{
    if(n == 3){
      results_binder <- binder(proxy_regulon, coexpression, operons, mu_phi=3, sigma_phi=0.01, is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
      df_binder <- data.frame("Regulator"=regulator, "TargetCandidate"=results_binder$target_candidate, "Score"=results_binder$theta_interval[, 2], "Interaction"=interaction_status, regulator_name=regulator_name, organism=organism, model="BINDER mu[phi]:=3, sigma[phi]:=0.01", ME=results_binder$ME, PE=results_binder$PE, CM=results_binder$CM, CP=results_binder$CP)
      all_df_binder <- rbind(all_df_binder, df_binder)
    }
    results_binder <- binder(proxy_regulon, coexpression, operons, is_coexpression=TRUE, chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
    df_binder <- data.frame("Regulator"=regulator, "TargetCandidate"=results_binder$target_candidate, "Score"=results_binder$theta_interval[, 2], "Interaction"=interaction_status, regulator_name=regulator_name, organism=organism, model="BINDER", ME=results_binder$ME, PE=results_binder$PE, CM=results_binder$CM, CP=results_binder$CP)
    all_df_binder <- rbind(all_df_binder, df_binder)

    results_non_auxiliary <- binder(proxy_regulon, coexpression, operons, is_coexpression=TRUE, model="non_auxiliary", chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
    df_non_auxiliary <- data.frame("Regulator"=regulator, "TargetCandidate"=results_non_auxiliary$target_candidate, "Score"=results_non_auxiliary$theta_interval[, 2], "Interaction"=interaction_status, regulator_name=regulator_name, organism=organism, model="Non-auxiliary", ME=results_non_auxiliary$ME, PE=results_non_auxiliary$PE, CM=results_non_auxiliary$CM, CP=results_non_auxiliary$CP)
    all_df_non_auxiliary <- rbind(all_df_non_auxiliary, df_non_auxiliary)

    results_deterministic <- binder(proxy_regulon, coexpression, operons, is_coexpression=TRUE, model="deterministic", chains=n_chains, iter=n_iter, warmup=n_warmup, seed=seed)
    df_deterministic <- data.frame("Regulator"=regulator, "TargetCandidate"=results_deterministic$target_candidate, "Score"=results_deterministic$theta_interval[, 2], "Interaction"=interaction_status, regulator_name=regulator_name, organism=organism, model="Deterministic", ME=results_deterministic$ME, PE=results_deterministic$PE, CM=results_deterministic$CM, CP=results_deterministic$CP)
    all_df_deterministic <- rbind(all_df_deterministic, df_deterministic)
  }
}

#########################
#########################

all_df <- rbind(all_df_binder[, 1:7], all_df_non_auxiliary[, 1:7], all_df_deterministic[, 1:7], all_df_genie)
all_df$combination <- paste(all_df$regulator_name, all_df$organism, all_df$model, sep="-")

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

base_palette <- rev(c("#2E0219", "#832232", "#CE8964", "#EFC69B", "#54b8d3", "grey"))

ggplot(all_df, aes(d=Interaction, m=Score, color=model, size=model)) +
  geom_abline(aes(intercept=0, slope=1), linetype=2) +
  geom_roc(n.cuts=0) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  facet_grid(regulator_name ~ organism) +
  theme_minimal() +
  theme(legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text=element_text(size=15)) +
  scale_color_manual(name="Model", labels=c("BINDER (all)", expression(paste("BINDER [Informative p(",phi,")]")),"BINDER","non-auxiliary","deterministic","GENIE3"), values=base_palette) +
  scale_size_manual(name="Model", labels=c("BINDER (all)", expression(paste("BINDER [Informative p(",phi,")]")),"BINDER","non-auxiliary","deterministic","GENIE3"), values=c(1,0.5,1,0.5,0.5,1)) +
  guides(shape=guide_legend(override.aes=list(size=1)))

ggplot(all_df) +
  geom_boxplot(aes(x=ifelse(Interaction == 1, "Interaction", "No Interaction"),y=Score,color=model,fill=model), alpha=0.7) +
  xlab("Interaction") +
  ylab(expression(paste(theta^0.5))) +
  facet_grid(regulator_name ~ organism) +
  theme_minimal() +
  theme(axis.text.y=element_text(angle=90), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text=element_text(size=15)) +
  scale_fill_manual(name="Model", labels=c("BINDER (all)", expression(paste("BINDER [Informative p(",phi,")]")),"BINDER","non-auxiliary","deterministic","GENIE3"), values=base_palette) +
  scale_color_manual(name="Model", labels=c("BINDER (all)", expression(paste("BINDER [Informative p(",phi,")]")),"BINDER","non-auxiliary","deterministic","GENIE3"), values=base_palette) +
  guides(shape=guide_legend(override.aes=list(size=1))) +
  coord_flip()

all_df <- rbind(all_df_binder, all_df_non_auxiliary, all_df_deterministic)
theta_score_difference <- data.frame(
  ME=ifelse(all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$ME == 1, 1, 0),
  PE=ifelse(all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$PE == 1, 1, 0),
  CM=all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$CM,
  CP=all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$CP,
  interaction=rep(ifelse(regulon_data$b_subtilis_lexA$interaction_status == 1, "Interaction", "No Interaction"), 3),
  score=all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$Score,
  model=all_df[all_df$regulator_name == "lexA" & all_df$organism == "B. subtilis" & all_df$model != "GENIE3", ]$model
)

ggplot(theta_score_difference[theta_score_difference$ME == 0 & theta_score_difference$PE == 0, ]) +
  geom_jitter(aes(x=CM, y=CP, fill=score, color=interaction), shape=22, size=3, width=0.1, height=0.1, alpha=0.75) +
  xlab("CM") +
  ylab("CP") +
  facet_grid(model ~ interaction) +
  theme_minimal() +
  theme(legend.title=element_text(size=12), legend.text=element_text(size=10), strip.text=element_text(size=15)) +
  scale_fill_gradientn(name=expression(paste(theta^0.5)), colours=c("#efefef","#7a7a7a"), limits=c(0,1)) +
  scale_color_manual(name="Validated Interaction", values=c("#832232","#EFC69B"))
