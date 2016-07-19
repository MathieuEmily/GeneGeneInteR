# Function that computes the SSI matrix between two given genes.
SSI.test <- function(Y, G1, G2){
  # Arguments checks
  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(G1) != "SnpMatrix" | class(G2) != "SnpMatrix") {
    stop("G1 and G2 arguments should be SnpMatrix objects.")
  } else if (nrow(G1) != nrow(G2)) {
    stop("G1 and G2 should have same rows count.")
  } else if (length(Y) != nrow(G1) | length(Y) != nrow(G2)) {
    stop("Response variable should be conformant with genes matrices rows number.")
  }

  snp.interaction <- matrix(0, nrow=ncol(G1), ncol=ncol(G2), dimnames=list(colnames(G1), colnames(G2)))
  SS.combn <- expand.grid(colnames(G1), colnames(G2))

  for (i in 1:nrow(SS.combn)) {
    S1.ind <- which(colnames(G1) == SS.combn[i, 1])
    S2.ind <- which(colnames(G2) == SS.combn[i, 2])
    snp.interaction[SS.combn[i, 1], SS.combn[i, 2]] <- SSI.BaseTest(Y, G1[, S1.ind], G2[, S2.ind])
  }

  return(snp.interaction)
}

# Function that computes interaction p-value between two
# SNPs. Called by SSI.test
SSI.BaseTest <- function(Y, S1, S2){
  # Arguments checks
  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(S1) != "SnpMatrix" | class(S2) != "SnpMatrix") {
    stop("S1 and S2 arguments should be SnpMatrix objects.")
  } else if (nrow(S1) != nrow(S2)) {
    stop("S1 and S2 should have same rows count.")
  } else if (length(Y) != nrow(S1) | length(Y) != nrow(S2)) {
    stop("Response variable should be conformant with snps matrices rows number.")
  }

  data <- data.frame(Y=Y, S1=as.numeric(S1), S2=as.numeric(S2))

  null.mod   <- glm(Y ~ S1 + S2, data=data, family="binomial")
  inter.mod <- glm(Y ~ S1*S2, data=data, family="binomial")

  comp.res <- anova(null.mod, inter.mod, test="Chisq")

  return(comp.res$'Pr(>Chi)'[2])
}
