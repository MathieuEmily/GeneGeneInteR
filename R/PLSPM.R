PLSPM.test <- function(Y, G1, G2, n.perm=500){

	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))
	

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(nlevels(as.factor(Y))!=2){
    stop("Y must be a factor with 2 levels, most likely 0 and 1.")
  } else if(!is(G1,"SnpMatrix")){
    stop("G1 must be a SnpMatrix object.")
  } else if(!is(G2,"SnpMatrix")){
    stop("G2 must be a SnpMatrix object")
  } else if(nrow(G1)!=nrow(G2)){
    stop("Both G1 and G2 must contain the same number of individuals.")
  } else if(length(Y)!=nrow(G1)){
    stop("Y and both SnpMatrix objects must contain the same number of individuals.")
  } else if (sum(is.na(G1))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(G2))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(Y))!=0) {
    stop("The response variable must be complete. No NAs are allowed.")
  }

  X1 <- as(G1, "numeric")
  X2 <- as(G2, "numeric")

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  X <- cbind(X1,X2)

  Gene1 <- c(0,0)
  Gene2 <- c(1,0)

  w1 <- which(Y==1)
  w0 <- which(Y==0)

  XCases <- X[w1,]
  XControls <- X[w0,]

  my.path <- rbind(Gene1,Gene2)
  my.blocks <- list(seq_len(ncol(X1)),(ncol(X1)+1):(ncol(X1)+ncol(X2)))
  my.modes = c("A", "A")

  mod1<-NULL;
  mod0<-NULL;

  try(mod1 <- plspm(XCases,my.path,my.blocks, modes = my.modes), silent=TRUE)
  if(is.null(mod1)){
  	warning("P-value could not be computed. NA returned")
	list.param <- list(n.perm=n.perm)
  	res <- list(statistic=NA,p.value=NA,method="Partial Least Squares Path Modeling",parameter=list.param)
	class(res) <- "GGItest"
	return(res)
	}

  try(mod0 <- plspm(XControls,my.path,my.blocks, modes = my.modes),silent=TRUE)
  if(is.null(mod0)){
  	warning("P-value could not be computed. NA returned")
	list.param <- list(n.perm=n.perm)
  	res <- list(statistic=NA,p.value=NA,method="Partial Least Squares Path Modeling",parameter=list.param)
	class(res) <- "GGItest"
	return(res)
	}

  beta1 <- mod1$inner_model[[1]][2,1]
  vbeta1 <- mod1$inner_model[[1]][2,2]^2
  beta0 <- mod0$inner_model[[1]][2,1]
  vbeta0 <- mod0$inner_model[[1]][2,2]^2

  U <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  
  U.perm <- rep(NA,times=n.perm)
  for (i in seq_len(n.perm)){
    restart<-TRUE
    while(restart){
	  	Y.perm <- sample(Y)
  		w1 <- which(Y.perm==1)
	  	w0 <- which(Y.perm==0)
  		XCases <- X[w1,]
	  	XControls <- X[w0,]
  		mod1<-NULL;
		mod0<-NULL;
		try(mod1 <- plspm(XCases,my.path,my.blocks, modes = my.modes), silent=TRUE)
#		if(is.null(mod1)){warning("P-value could not be computed. NA returned");return(NA)}
		try(mod0 <- plspm(XControls,my.path,my.blocks, modes = my.modes),silent=TRUE)
#		if(is.null(mod0)){warning("P-value could not be computed. NA returned");return(NA)}
		if (!is.null(mod1) & !is.null(mod0)){restart <- FALSE}
	}
	beta1 <- mod1$inner_model[[1]][2,1]
	vbeta1 <- mod1$inner_model[[1]][2,2]^2
	beta0 <- mod0$inner_model[[1]][2,1]
	vbeta0 <- mod0$inner_model[[1]][2,2]^2
	
	U.perm[i] <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  }
  
  #pval <- 2*(1-pnorm(abs(U)))
	  pval <- mean(abs(U.perm) > abs(U))
	  stat <- U
	names(stat)="U"
#	list.param <- list(n.perm=n.perm)
#	res <- list(statistic=stat,p.value=pval,method="Partial Least Squares Path Modeling",parameter=list.param)
#	class(res) <- "GGItest"
 # return(res)
  
    null.value <- 0
names(null.value) <- "U"
estimate <- c(beta0, beta1)
names(estimate) <- c("beta0","beta1")
parameters <- n.perm
names(parameters) <- "n.perm"
	res <- list(
		null.value=null.value,
		alternative="two.sided",
		method="Gene-based interaction based on Partial Least Squares Path Modeling",
		estimate= estimate,
		data.name=paste(Y.arg," and  (",G1.arg," , ",G2.arg,")",sep=""),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)


#  return(pval)
}


plspm <-
  function(Data, path_matrix, blocks, modes = NULL, scaling = NULL,  
           scheme = "centroid", scaled = TRUE, tol = 0.000001, maxiter = 100, 
           plscomp = NULL, boot.val = FALSE, br = NULL, 
           dataset = TRUE)
  {
    # =======================================================================
    # checking arguments
    # =======================================================================
    valid = check_args(Data=Data, path_matrix=path_matrix, blocks=blocks, 
                       scaling=scaling, modes=modes, scheme=scheme, 
                       scaled=scaled, tol=tol, maxiter=maxiter, 
                       plscomp=plscomp, boot.val=boot.val, br=br, 
                       dataset=dataset)
    
    Data = valid$Data
    path_matrix = valid$path_matrix
    blocks = valid$blocks
    specs = valid$specs
    boot.val = valid$boot.val
    br = valid$br
    dataset = valid$dataset
    
    # =======================================================================
    # Preparing data and blocks indexification
    # =======================================================================
    # building data matrix 'MV'
    MV = get_manifests(Data, blocks)
    check_MV = test_manifest_scaling(MV, specs$scaling)
    # generals about obs, mvs, lvs
    gens = get_generals(MV, path_matrix)
    # blocks indexing
    names(blocks) = gens$lvs_names
    block_sizes = lengths(blocks)
    blockinds = indexify(blocks)
    
    # transform to numeric if there are factors in MV
    if (test_factors(MV)) {
      numeric_levels = get_numerics(MV)
      MV = numeric_levels$MV
      categories = numeric_levels$categories
    }  
    # apply corresponding treatment (centering, reducing, ranking)
    X = get_treated_data(MV, specs)
    
    # =======================================================================
    # Outer weights and LV scores
    # =======================================================================
    metric = get_metric(specs$scaling)
    if (metric) {
      # object 'weights' contains outer w's, W, ODM, iter
      weights = get_weights(X, path_matrix, blocks, specs)
      ok_weights = test_null_weights(weights, specs)
      outer_weights = weights$w
      LV = get_scores(X, weights$W)
    } else {
      # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
      weights = get_weights_nonmetric(X, path_matrix, blocks, specs)
      ok_weights = test_null_weights(weights, specs)
      outer_weights = weights$w
      LV = weights$Y
      X = weights$QQ  # quantified MVs
      colnames(X) = gens$mvs_names
    }
    
    # =======================================================================
    # Path coefficients and total effects
    # =======================================================================
    inner_results = get_paths(path_matrix, LV)
    inner_model = inner_results[[1]]
    Path = inner_results[[2]]
    R2 = inner_results[[3]]
    Path_effects = get_effects(Path)
    
    # =======================================================================
    # Outer model: loadings, communalities, redundancy, crossloadings
    # =======================================================================
    xloads = cor(X, LV, use = 'pairwise.complete.obs')
    loadings = rowSums(xloads * weights$ODM)
    communality = loadings^2
    R2_aux = rowSums(weights$ODM %*% diag(R2, gens$lvs, gens$lvs))
    redundancy = communality * R2_aux
    crossloadings = data.frame(xloads, row.names=seq_len(gens$mvs))
    crossloadings$name = factor(gens$mvs_names, levels = unique(gens$mvs_names))
    crossloadings$block = factor(rep(gens$lvs_names, block_sizes),
                                 levels = gens$lvs_names)
    crossloadings = crossloadings[,c('name','block',colnames(xloads))]
    
    # outer model data frame
    outer_model = data.frame(
      name = factor(gens$mvs_names, levels = unique(gens$mvs_names)),
      block = factor(rep(gens$lvs_names, block_sizes),
                     levels = gens$lvs_names),
      weight = outer_weights, 
      loading = loadings, 
      communality = communality,
      redundancy = redundancy,
      #row.names = 1:gens$mvs)
      row.names = seq_len(gens$mvs))
    
    # Unidimensionality
    unidim = get_unidim(DM = MV, blocks = blocks, modes = specs$modes)
    
    # Summary Inner model
    inner_summary = get_inner_summary(path_matrix, blocks, specs$modes,
                                      communality, redundancy, R2)
    
    # GoF Index
    gof = get_gof(communality, R2, blocks, path_matrix)
    
    # =======================================================================
    # Results
    # =======================================================================
    # deliver dataset?
    if (dataset) data = MV else data = NULL
    # deliver bootstrap validation results? 
    bootstrap = FALSE
    if (boot.val) 
    {
      if (nrow(X) <= 10) {
        warning("Bootstrapping stopped: very few cases.") 
      } else { 
        bootstrap = get_boots(MV, path_matrix, blocks, specs, br)
      }
    }
    
    # list with model specifications
    model = list(IDM=path_matrix, blocks=blocks, specs=specs,
                 iter=weights$iter, boot.val=boot.val, br=br, gens=gens)
    
    ## output
    res = list(outer_model = outer_model, 
               inner_model = inner_model,
               path_coefs = Path, 
               scores = LV,
               crossloadings = crossloadings, 
               inner_summary = inner_summary, 
               effects = Path_effects,
               unidim = unidim, 
               gof = gof, 
               boot = bootstrap, 
               data = data,
               manifests = X,
               model = model)
    class(res) = "plspm"
    return(res)
  }


#'@S3method print plspm
print.plspm <- function(x, ...)
{
  cat("Partial Least Squares Path Modeling (PLS-PM)", "\n")
  cat(rep("-", 45), sep="")
  cat("\n   NAME            ", "DESCRIPTION")  
  cat("\n1  $outer_model    ", "outer model")
  cat("\n2  $inner_model    ", "inner model")
  cat("\n3  $path_coefs     ", "path coefficients matrix")
  cat("\n4  $scores         ", "latent variable scores")
  if (!inherits(x, "plspm.fit"))
  {
    cat("\n5  $crossloadings  ", "cross-loadings")
    cat("\n6  $inner_summary  ", "summary inner model")
    cat("\n7  $effects        ", "total effects")
    cat("\n8  $unidim         ", "unidimensionality")
    cat("\n9  $gof            ", "goodness-of-fit")
    cat("\n10 $boot           ", "bootstrap results")
    cat("\n11 $data           ", "data matrix")
  }
  cat("\n")
  cat(rep("-", 45), sep="")
  cat("\nYou can also use the function 'summary'", "\n\n")    
  invisible(x)
}


check_args <- 
  function(Data, path_matrix, blocks, scaling, modes, scheme,
           scaled, tol, maxiter, plscomp, boot.val, br, dataset)
  {
    # check definitions
    Data = check_data(Data)
    path_matrix = check_path(path_matrix)
    blocks = check_blocks(blocks, Data)
    specs = check_specs(blocks, scaling, modes, scheme, scaled, 
                        tol, maxiter, plscomp)
    boot_args = check_boot(boot.val, br)
    if (!is.logical(dataset)) dataset = TRUE
    
    # check congruence between inner model and outer model
    good_model = check_model(path_matrix, blocks)
    
    # list with verified arguments
    list(Data = Data,
         path_matrix = path_matrix,
         blocks = blocks,
         specs = specs,
         boot.val = boot_args$boot.val,
         br = boot_args$br, 
         dataset = dataset)
  }

check_data <- function(Data)
{
  if (is_not_tabular(Data))
    stop("\nInvalid 'Data'. Must be a matrix or data frame.")
  
  if (is.matrix(Data) && !is.numeric(Data))
    stop("\nInvalid 'Data' matrix. Must be a numeric matrix.")
  
  if (nrow(Data) == 1)
    stop("\nCannot work with only one row in 'Data'")
  
  if (ncol(Data) == 1)
    stop("\nCannot work with only one column in 'Data'")
  
  if (lacks_rownames(Data))
    rownames(Data) = seq_len(nrow(Data))
  
  if (lacks_colnames(Data)) 
    colnames(Data) = paste("MV", seq_len(ncol(Data)), sep="")
  
  # return
  Data
}

check_path <- function(path_matrix)
{
  if (is_not_matrix(path_matrix))
    stop("\n'path_matrix' must be a matrix.")
  
  if (!is_square_matrix(path_matrix))
    stop("\n'path_matrix' must be a square matrix.")
  
  if (nrow(path_matrix) == 1)
    stop("\n'path_matrix' must have more than one row")
  
  if (!is_lower_triangular(path_matrix))
    stop("\n'path_matrix' must be a lower triangular matrix")
  
  
  for (j in seq_len(ncol(path_matrix))) 
  {
    for (i in seq_len(nrow(path_matrix)))
    {
      if (length(intersect(path_matrix[i,j], c(1,0))) == 0)
        stop("\nElements in 'path_matrix' must be '1' or '0'")
    }      
  }
  
  if (lacks_dimnames(path_matrix)) {
    LV_names = paste("LV", seq_len(ncol(path_matrix)), sep = "")
    dimnames(path_matrix) = list(LV_names, LV_names)
  }
  if (has_rownames(path_matrix) && lacks_colnames(path_matrix)) {
    colnames(path_matrix) = rownames(path_matrix)
  }
  if (has_colnames(path_matrix) && lacks_rownames(path_matrix)) {
    rownames(path_matrix) = colnames(path_matrix)
  }
  
  # return
  path_matrix
}

check_blocks <- function(blocks, Data)
{
  if (!is.list(blocks))
    stop("\n'blocks' must be a list.")
  
  # no duplicated elements within each block
  mvs_duplicated = unlist(lapply(blocks, duplicated))
  if (any(mvs_duplicated))
    stop("\nWrong 'blocks'. Duplicated variables in a block are not allowed")
  
  # all elements in blocks of same mode
  mvs_mode = unique(unlist(lapply(blocks, mode)))
  if (length(mvs_mode) > 1)
    stop("\nAll elements in 'blocks' must have the same mode")
  
  # check indices inside columns range of Data
  if (mvs_mode == "numeric") {
    blocks_in_data = match(unlist(blocks), seq_len(ncol(Data)))
    if (any(is.na(blocks_in_data)))
      stop("\nIndices in 'blocks' outside the number of columns in 'Data'")
  }
  
  # convert character blocks to numeric blocks
  if (mvs_mode == "character") {
    data_names = colnames(Data)
    matched_names = match(unlist(blocks), data_names)
    if (any(is.na(matched_names))) {
      bad_names = unlist(blocks)[is.na(matched_names)]
      stop(sprintf("\nUnrecognized name in 'blocks': '%s'", bad_names))        
    }
    blocks = lapply(blocks, function(x, y) match(x, y), data_names)
  }
  
  # output
  blocks
}

check_boot <- function(boot.val, br)
{
  if (!is.logical(boot.val)) boot.val = FALSE
  
  if (boot.val) {
    if (!is.null(br)) {
      if(!is_positive_integer(br) || length(br) != 1L || br < 10) {
        warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
        br = 100
      } 
    } else
      br = 100
  }
  
  # return
  list(boot.val = boot.val, br = br)
}


check_model <- function(path_matrix, blocks)
{    
  # compatibility between path_matrix and blocks
  if (length(blocks) != nrow(path_matrix))
    stop("\nNumber of rows in 'path_matrix' different from length of 'blocks'.")
  
  # output
  TRUE
}

is_tabular <- function(x) {
  if (is.matrix(x) | is.data.frame(x)) {
    TRUE
  } else FALSE  
}

is_numeric_tabular <- function(x) {
  if (is_numeric_matrix(x) | is_numeric_dataframe(x)) {
    TRUE
  } else FALSE  
}

is_string_tabular <- function(x) {
  if (is_string_matrix(x) | is_string_dataframe(x)) {
    TRUE
  } else FALSE  
}

is_not_tabular <- function(x) {
  !is_tabular(x)  
}

has_rownames <- function(x) {
  if (!is.null(rownames(x))) TRUE else FALSE
}

has_colnames <- function(x) {
  if (!is.null(colnames(x))) TRUE else FALSE
}

has_dimnames <- function(x) {
  if (!is.null(dimnames(x))) TRUE else FALSE
}

lacks_rownames <- function(x) {
  !has_rownames(x)
}

lacks_colnames <- function(x) {
  !has_colnames(x)
}

lacks_dimnames <- function(x) {
  !has_dimnames(x)
}

is_matrix <- function(x) {
  is.matrix(x)
}

is_numeric_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.numeric(x)
}

is_string_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.character(x)
}

is_logical_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.logical(x)
}

is_not_matrix <- function(x) {
  !is_matrix(x)
}

is_square_matrix <- function(x) {
  if (is.matrix(x)) {
    if (nrow(x) == ncol(x)) TRUE else FALSE      
  } else FALSE
}

is_not_square_matrix <- function(x) {
  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) TRUE else FALSE      
  } else TRUE
}

is_square_numeric_matrix <- function(x) {
  if (is_numeric_matrix(x)) {
    if (nrow(x) == ncol(x)) TRUE else FALSE      
  } else FALSE
}

is_not_square_numeric_matrix <- function(x) {
  !is_square_numeric_matrix(x)
}

is_diagonal <- function(x) {
  if (is_square_matrix(x)) {
    above = sum(x[upper.tri(x)])
    below = sum(x[lower.tri(x)])
    if (above > 0 | below > 0) FALSE else TRUE      
  } else FALSE
}

is_not_diagonal <- function(x) {
  !is_diagonal(x)
}

is_lower_triangular <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    all(x[upper.tri(x, diag = diag)] == 0)
  } else FALSE
}

is_upper_triangular <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    all(x[lower.tri(x, diag = diag)] == 0)
  } else FALSE
}

is_triangular_matrix <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    is_lower_triangular(x) | is_upper_triangular(x)
  } else FALSE
}

check_specs <-
  function(blocks, scaling, modes, scheme, scaled, tol, maxiter, plscomp)
  {
    check_scale = check_scaling(scaling, scaled, blocks)
    
    # output
    list(scaling = check_scale$scaling,
         modes = check_modes(modes, blocks),
         scheme = check_scheme(scheme),
         scaled = check_scale$scaled,
         tol = check_tol(tol),
         maxiter = check_maxiter(maxiter),
         plscomp = check_plscomp(plscomp, scaling, modes))  
  }

check_scaling <- function(scaling, scaled, blocks)
{
  # if scaling is present
  if (!is.null(scaling)) 
  {
    # make sure scaling is a list
    if (!is.list(scaling))
      stop("\nInvalid 'scaling'. Must be a list.")
    # compatibility between blocks and scaling
    if (length(blocks) != length(scaling)) {
      stop("\nLength of 'scaling' differs from length of 'blocks'.")
    }
    if (!identical(lengths(blocks), lengths(scaling))) {
      stop("\nLengths of 'scaling' differs from lengths of 'blocks'.")
    }
    
    # string manipulation of elements in 'scaling'
    scaling_aux = tolower(unlist(scaling))
    scaling_aux = substr(scaling_aux, start=1, stop=3)
    
    # are there any unrecognized scaling types?
    bad_scale <- !(scaling_aux %in% c("num", "raw", "ord", "nom"))
    if (any(bad_scale)) {
      bad = unlist(scaling)[bad_scale]
      stop(sprintf("\nSorry. Unrecognized scaling type: '%s'", bad))
    }
    
    # set all numeric when mixing only 'num' and 'raw'
    if (!any(scaling_aux %in% c('ord', 'nom'))) {
      num_and_raw = intersect(scaling_aux, c("num", "raw"))
      if (length(num_and_raw) > 1) {
        scaling = lapply(blocks, function(x) rep("num", length(x)))
      }
      if (all(unlist(scaling) == "num")) scaled = TRUE
      if (all(unlist(scaling) == "raw")) scaled = FALSE    
    }
    # final scaling
    scaling = lapply(scaling, function(x) substr(tolower(x), start=1, stop=3))
  } else {
    if (!is.logical(scaled)) scaled = TRUE
  }
  
  # output
  list(scaling = scaling, scaled = scaled)
}

check_modes <- function(modes, blocks)
{
  # default modes    
  if (is.null(modes)) {
    modes = rep("A", length(blocks))
  } 
  
  if (length(blocks) != length(modes)) {
    warning("Warning: length of 'modes' different from length of 'blocks'")
    message("Default modes 'A' is used")
    modes = rep("A", length(blocks))
  }
  
  # are there any unrecognized modes?
  modes = toupper(modes)
  bad_modes <- !(modes %in% c("A", "B", "NEWA", "PLSCOW", "PLSCORE"))
  if (any(bad_modes)) {
    bad = modes[bad_modes]
    stop(sprintf("\nSorry. Unrecognized mode: '%s'", bad))
  }
  
  # cannot mix modes "A" and "newA"
  mixed_modes = intersect(modes, c("A", "NEWA"))
  if (length(mixed_modes) > 1) {
    stop("\nSorry. Can't work with both modes 'A' and 'newA'")
  }  
  
  # output
  modes
}

check_scheme <- function(scheme)
{
  # some string manipulations
  if (!is.character(scheme)) scheme = as.character(scheme)
  scheme = tolower(scheme)
  
  SCHEMES = c("centroid", "factorial", "path")
  scheme_match = pmatch(scheme, SCHEMES)
  if (is.na(scheme_match)) {
    warning("\nInvalid 'scheme'. Default 'scheme=centroid' is used.")   
    scheme = "centroid"
  } else {
    scheme = SCHEMES[scheme_match]
  }
  
  # output
  scheme
}

check_tol <- function(tol)
{
  if (mode(tol) != "numeric" || length(tol) != 1L || 
      tol <= 0 || tol > 0.001) {
    warning("Warning: Invalid 'tol'. Default 'tol=0.000001' is used.")   
    tol = 0.000001
  } 
  
  # output
  tol
}

check_maxiter <- function(maxiter)
{
  if (!is_positive_integer(maxiter) || maxiter < 100) {
    warning("Warning: Invalid 'maxiter'. Default 'maxiter=100' is used.")   
    maxiter = 100
  } 
  
  # output
  maxiter
}

check_plscomp <- function(plscomp, scaling, modes)
{
  if (is.null(scaling)) plscomp = NULL
  
  if (!is.null(scaling)) {
    if (is.null(plscomp)) {
      plscomp = rep(1, length(scaling))      
    } else {
      if (any(!is_positive_integer(plscomp)))
        stop("\n'plscomp' must be an integer vector")
      if (length(scaling) != length(plscomp))
        stop("\nlength of 'plscomp' differs from number of blocks")
      
      plscomp = as.integer(plscomp)
      for (j in seq_len(length(plscomp)))
      {
        # plscomp's cannot exceed number of variables in block j
        if (plscomp[j] > length(scaling[[j]]))
        {
          stop(sprintf("%s %d %s", "element", j, 
                       "in 'plscomp' exceeds number of variables"))          
        }
        # make sure plscomp[j]=1 when mode "NEWA"
        if (modes[j] == "NEWA") plscomp[j] = 1
      }
    }
  }
  
  # output
  plscomp
}

is_integer <- function(x) {
  UseMethod("is_integer", x)
}

#' @S3method is_integer default
is_integer.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_integer factor
is_integer.factor <- function(x) {
  FALSE
}

#' @S3method is_integer numeric
is_integer.numeric <- function(x) {
  (x %% 1) == 0
}

is_not_integer <- function(x) {
  !is_integer(x)
}

is_positive_integer <- function(x) {
  (is_positive(x) & is_integer(x))
}

is_negative_integer <- function(x) {
  (is_negative(x) & is_integer(x))
}

is_positive <- function(x) {
  UseMethod("is_positive", x)
}

#' @S3method is_positive default
is_positive.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_positive numeric
is_positive.numeric <- function(x) {
  x > 0
}

#' @S3method is_positive matrix
is_positive.matrix <- function(x) {
  x > 0
}

#' @S3method is_positive factor
is_positive.factor <- function(x) {
  FALSE
}

is_not_positive <- function(x) {
  !is_positive(x)
}

is_negative <- function(x) {
  UseMethod("is_negative", x)
}

#' @S3method is_negative default
is_negative.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_negative numeric
is_negative.numeric <- function(x) {
  x < 0
}

#' @S3method is_negative matrix
is_negative.matrix <- function(x) {
  x < 0
}

#' @S3method is_negative factor
is_negative.factor <- function(x) {
  FALSE
}

is_not_negative <- function(x) {
  !is_negative(x)
}

get_manifests <- function(Data, blocks)
{
  # building data matrix 'MV'
  MV = Data[,unlist(blocks)]
  
  # add row and column names
  mvs_names = colnames(Data)[unlist(blocks)]
  dimnames(MV) = list(rownames(Data), mvs_names)
  
  # output
  MV
}

test_manifest_scaling <- function(MV, scaling)
{
  # to be used when MV is a data.frame containing factors
  if (is.data.frame(MV)) {
    mvs_class = vapply(MV, class)
    mvs_as_factors <- mvs_class == "factor"
    # if there are MVs as factors, check right scaling
    if (sum(mvs_as_factors) > 0) {
      factors_scaling = unlist(scaling)[mvs_as_factors]
      
      # factors can't be numeric or raw
      if (any(factors_scaling %in% c("num", "raw")))
        stop("\n'Data' contains factors that can't have metric scaling")
      
      # unordered factors must be nominal
      if (sum(factors_scaling == "ord") == 1) {
        if (!is.ordered(MV[,mvs_as_factors]))
          stop("\nUnordered factors in 'Data' can't have ordinal scaling")
      } 
      if (sum(factors_scaling == "ord") > 1)  {
        unordered = !apply(MV[,mvs_as_factors], 2, is.ordered)
        if (sum(unordered) > 0)
          stop("\nUnordered factors in 'Data' can't have ordinal scaling")
      }    
    }    
  }
  TRUE
}

get_generals <- function(MV, path_matrix)
{
  list(obs = nrow(MV),
       obs_names = rownames(MV),
       mvs = ncol(MV),
       mvs_names = colnames(MV),
       lvs = nrow(path_matrix),
       lvs_names = rownames(path_matrix))
}

indexify <- function(x, out) {
  UseMethod("indexify", x)  
}


#' @S3method indexify default
indexify.default <- function(x, ...)
{
  if (!is_numeric_vector(x) || !list_of_vectors(x))
    stop("\n'indexify()' requires a numeric vector or a list of vectors")
}


#' @S3method indexify numeric
indexify.numeric <- function(x, out = "vector")
{
  if (!is_numeric_vector(x))
    stop("\n'indexify()' requires a numeric vector")
  
  if (out == "vector")
    rep(seq_along(x), x)
  else mapply(rep, seq_along(x), x)
}


#' @S3method indexify list
indexify.list <- function(x, out = "vector")
{
  if (!list_of_vectors(x))
    stop("\n'indexify()' requires a list of vectors")
  
  aux = unlist(lapply(x, length))
  if (out == "vector")
    rep(seq_along(aux), aux)
  else mapply(rep, seq_along(aux), aux)
}

list_of_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is.vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_numeric_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_numeric_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_string_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_string_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_logical_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_logical_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

test_factors <- function(DF)
{
  contains_factors = FALSE
  # to be used when DF is a data.frame containing factors
  if (is.data.frame(DF)) {
    mvs_class = vapply(DF, class)
    mvs_as_factors <- mvs_class == "factor"
    # tell me if there are MVs as factors
    if (sum(mvs_as_factors) > 0)
      contains_factors = TRUE
  }
  contains_factors
}

get_numerics <- function(MV)
{  
  mvs_class = vapply(MV, class)
  mvs_as_factors <- mvs_class == "factor"
  categorical = which(mvs_as_factors)
  categories = vector(mode="list", ncol(MV))
  
  # only keep levels of categorical mvs in 'factor' format
  for (j in seq_along(categorical)) {
    # keep original levels
    categories[[categorical[j]]] = levels(MV[,categorical[j]])
    # convert to numeric
    MV[,categorical[j]] = as.numeric(MV[,categorical[j]])
  }    
  
  # output
  list(MV=MV, categories=categories)
}

get_treated_data <- function(MV, specs)
{
  metric = get_metric(specs$scaling)
  
  if (metric) {
    # standardize if all numeric
    if (specs$scaled) {
      sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
      X = scale(MV, scale=sd.X)
    } else {
      # center if all raw
      X = scale(MV, scale=FALSE)
    }
  } else {
    # standardize
    #    sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
    #    X = scale(MV, scale=sd.X)
    X = scale(MV) / sqrt((nrow(MV)-1)/nrow(MV))
    # all variables quantified as "ord"/"nom" are pre-treated with get_rank
    scaling = unlist(specs$scaling)
    nominal_ordinal = which(scaling %in% c("ord", "nom"))
    
    for (j in nominal_ordinal) {
      X[,j] = get_rank(X[,j])
    }
  }
  
  # add names
  dimnames(X) = dimnames(MV)
  # output
  X
}


get_metric <- function(scaling)
{
  # metric case
  if (is.null(scaling)) metric = TRUE else metric = FALSE
  
  # output
  metric
}

get_weights <- function(X, path_matrix, blocks, specs)
{
  lvs = nrow(path_matrix)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)-1) / nrow(X))   # std.dev factor correction
  blockinds = indexify(blocks)
  
  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  W = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w_old = rowSums(W)    
  iter = 1
  
  repeat 
  {
    # external estimation of LVs 'Y'
    Y = X %*% W
    Y = scale(Y) * sdv
    # matrix of inner weights 'e' 
    E <- switch(specs$scheme, 
                "centroid" = sign(cor(Y) * (path_matrix + t(path_matrix))),
                "factorial" = cor(Y) * (path_matrix + t(path_matrix)),
                "path" = get_path_scheme(path_matrix, Y))
    # internal estimation of LVs 'Z'
    Z = Y %*% E  
    #    Z = Z %*% diag(1/(apply(Z,2,sd)*sdv), lvs, lvs)  # scaling Z
    # computing outer weights 'w'
    for (j in seq_len(lvs))
    {
      if (specs$modes[j] == "A")
        W[blockinds==j,j] <- (1/nrow(X)) * Z[,j] %*% X[,blockinds==j] 
      if (specs$modes[j] == "B")
        W[blockinds==j,j] <- solve.qr(qr(X[,blockinds==j]), Z[,j])
    }
    w_new = rowSums(W)                
    w_dif = sum((abs(w_old) - abs(w_new))^2) 
    if (w_dif < specs$tol || iter == specs$maxiter) break
    w_old = w_new
    iter = iter + 1
  } # end repeat       
  
  # preparing results
  if (iter == specs$maxiter) {
    results = NULL
  } else {
    W = W %*% diag(1/(apply(X %*% W, 2, sd)*sdv), lvs, lvs)
    w_new = rowSums(W)                
    names(w_new) = colnames(X)
    dimnames(W) = list(colnames(X), rownames(path_matrix))    
    results = list(w = w_new, W = W, ODM = ODM, iter = iter)
  }
  # output
  return(results)
}


vector_to_dummy <- function(avector)
{
  if (!is_numeric_vector(avector))
    stop("\n'vector_to_dummy()' requires a numeric vector")
  
  num_rows = sum(avector)
  num_cols = length(avector)
  
  # starting-and-ending positions
  start_end = from_to(avector)
  from = start_end$from
  to = start_end$to
  
  # build dummy matrix
  dummy_matrix = matrix(0, num_rows, num_cols)
  for (k in seq_len(num_cols)) {
    dummy_matrix[from[k]:to[k],k] = 1
  }
  dummy_matrix
}

list_to_dummy <- function(alist)
{
  if (!list_of_vectors(alist))
    stop("\n'list_to_dummy()' requires a list of vectors")
  
  aux = lengths(alist)
  to = cumsum(aux)
  from = to - aux + 1
  dummy_matrix = matrix(0, sum(aux), length(aux))
  for (j in seq_along(aux))
    dummy_matrix[from[j]:to[j], j] = 1
  dummy_matrix
}

test_null_weights <- function(wgs, specs) 
{
  if (is.null(wgs)) {
    print(paste("Iterative process is non-convergent with 'maxiter'=", 
                specs$maxiter, " and 'tol'=", specs$tol, sep=""))
    message("Algorithm stops") 
    stop("")
  }
  TRUE
}

get_scores <- function(X, W) 
{
  lvs = ncol(W)
  # correlations between MVs and LVs
  LV = X %*% W
  cor.XY = cor(X, LV)
  # sign ambiguity
  ODM = W
  ODM[W != 0] = 1
  w_sign = sign(colSums(sign((cor.XY * ODM))))
  if (any(w_sign <= 0)) {
    w_sign[w_sign == 0] = -1
    # scores
    LV = LV %*% diag(w_sign, lvs, lvs)    
  }
  colnames(LV) = colnames(W)
  LV
}

get_paths <-  function(path_matrix, Y_lvs, full=TRUE)
{
  lvs_names = colnames(path_matrix)
  endogenous = as.logical(rowSums(path_matrix))
  num_endo = sum(endogenous)
  results = as.list(seq_len(num_endo))
  Path = path_matrix
  residuals = as.list(seq_len(num_endo))
  R2 = rep(0, nrow(path_matrix))
  
  for (aux in seq_len(num_endo))
  {
    # index for endo LV
    k1 <- which(endogenous)[aux]
    # index for indep LVs
    k2 = which(path_matrix[k1,] == 1)
    
    path_lm = summary(lm(Y_lvs[,k1] ~ Y_lvs[,k2]))
    Path[k1,k2] = path_lm$coef[-1,1]
    residuals[[aux]] = path_lm$residuals  
    R2[k1] = path_lm$r.squared
    inn_val = c(path_lm$r.squared, path_lm$coef[,1])
    # ----- NEW results
    inn_labels = c("Intercept", names(k2))
    rownames(path_lm$coefficients) = inn_labels
    results[[aux]] <- path_lm$coefficients
    # ----- OLD results
    # inn_lab = c("R2", "Intercept", 
    # paste(rep("path_",length(k2)),names(k2),sep=""))
    # names(inn_val) = NULL
    # results[[aux]] <- data.frame(concept=inn_lab, value=round(inn_val, 4))
  }
  names(results) = lvs_names[endogenous]  
  names(R2) = lvs_names
  
  # output
  list(results, Path, R2, residuals)
}

get_effects <- function(Path)
{
  # how many latent variables and their names
  lvs = nrow(Path)
  lvs_names = rownames(Path)
  
  # list for storing effects
  path_effects = as.list(seq_len(lvs-1))
  path_effects[[1]] = Path
  
  # when only 2 lvs
  if (lvs == 2) {
    indirect_paths = matrix(c(0,0,0,0), 2, 2)
    total_paths = Path
  } else {
    # when more than 2 lvs
    for (k in 2:(lvs-1)) {
      path_effects[[k]] = path_effects[[k-1]] %*% Path        
    }
    indirect_paths = matrix(0, lvs, lvs)
    for (k in 2:length(path_effects)) {
      indirect_paths = indirect_paths + path_effects[[k]]        
    }
    total_paths = Path + indirect_paths
  }
  
  # initialize
  efs_labels <- direct <- indirect <- total <- NULL
  for (j in seq_len(lvs)) {
    for (i in j:lvs) {
      if (i != j) {
        efs_labels = c(efs_labels, paste(lvs_names[j], "->", lvs_names[i]))
        direct = c(direct, Path[i,j])
        indirect = c(indirect, indirect_paths[i,j])
        total = c(total, total_paths[i,j])
      }
    }
  }
  
  # output
  data.frame(relationships = efs_labels, 
             direct = direct, 
             indirect = indirect, 
             total = total)
}

get_unidim <- function(DM, blocks, modes)
{
  # inputs setting
  lvs = length(blocks) 
  lvs_names = names(blocks)
  blockinds = indexify(blocks)
  block_sizes = lengths(blocks)
  #  blocklist = unlist(lapply(block_sizes, function(x) rep(x, x)))
  obs = nrow(DM)
  sdvf = sqrt((nrow(DM)-1) / nrow(DM)) 
  # missing data flags
  missing_data = vapply(DM, is_missing)
  
  # Unidimensionality
  Alpha = rep(1, lvs)    # Cronbach's Alpha for each block
  Rho = rep(1, lvs)      # D.G. Rho for each block
  eig.1st = rep(1, lvs)  # first eigenvalue
  eig.2nd = rep(0, lvs)  # second eigenvalue
  
  # calculate indices
  for (aux in seq_len(lvs))
  {
    if (any(missing_data[blockinds == aux]))
    {
      Alpha[aux] = NA
      Rho[aux] = NA
      eig.1st[aux] = NA
      eig.2nd[aux] = NA
    } else {
      if (block_sizes[aux] != 1) 
      { 
        # scaling data
        DM.block = DM[,blockinds==aux]
        stdev.X = apply(DM.block, 2, sd) * sdvf 
        X_uni = scale(DM.block, scale=stdev.X)
        if (nrow(X_uni) < ncol(X_uni)) {   # more columns than rows
          acp = princomp(t(X_uni)) 
          X.rho = t(X_uni)
        } else {   # more rows than columns
          acp = princomp(X_uni)
          X.rho = X_uni
        }
        if (modes[aux] == "A") 
        {
          p = ncol(X_uni)
          # cronbach's alpha
          a.denom = var(rowSums(X_uni)) * sdvf^2
          a.numer = 2 * sum(cor(X_uni)[lower.tri(cor(X_uni))])
          alpha = (a.numer / a.denom) * (p / (p - 1))
          Alpha[aux] <- ifelse(alpha < 0, 0, alpha)
          # dillon-goldstein rho
          numer_rho <- colSums(cor(X.rho, acp$scores[,1]))^2
          denom_rho <- numer_rho + (p - colSums(cor(X.rho, acp$scores[,1])^2) )
          Rho[aux] <- numer_rho / denom_rho
        } else {  # modes[aux]=="B"
          Alpha[aux] = 0
          Rho[aux] = 0
        }
        eig.1st[aux] = acp$sdev[1]^2
        eig.2nd[aux] = acp$sdev[2]^2
      }      
    }
  }
  unidim = data.frame(Mode = modes, 
                      MVs = block_sizes,
                      C.alpha = Alpha, 
                      DG.rho = Rho,
                      eig.1st, 
                      eig.2nd)
  rownames(unidim) = lvs_names
  return(unidim)
}

is_missing <- function(x) 
{
  if (is.matrix(x)) {
    num_miss = apply(x, 2, function(u) sum(is.na(u)))    
  } else {
    num_miss = sum(is.na(x))
  }
  as.logical(sum(num_miss))
}

get_inner_summary <- 
  function(path_matrix, blocks, modes, communality, redundancy, R2)
  {
    blocklist = indexify(blocks)  
    exo_endo = rep("Exogenous", nrow(path_matrix))
    exo_endo[rowSums(path_matrix) != 0] = "Endogenous"
    avg_comu = rep(0, nrow(path_matrix))
    avg_redu = rep(0, nrow(path_matrix))
    AVE = rep(0, nrow(path_matrix))
    
    for (k in seq_len(nrow(path_matrix)))
    {
      avg_comu[k] = mean(communality[blocklist == k])
      avg_redu[k] = mean(redundancy[blocklist == k])
      if (modes[k] == "A")
      {
        ave_num = sum(communality[blocklist == k])
        ave_denom = ave_num + sum(1 - (communality[blocklist == k]))
        AVE[k] = ave_num / ave_denom
      }
    }
    
    # output
    data.frame(Type = exo_endo, 
               R2 = R2, 
               Block_Communality = avg_comu, 
               Mean_Redundancy = avg_redu, 
               AVE = AVE,
               row.names = rownames(path_matrix))
  }

get_gof <- function(comu, R2, blocks, path_matrix)
{
  lvs = nrow(path_matrix)
  blocklist = indexify(blocks)  
  endo = rowSums(path_matrix)
  endo[endo != 0] = 1  
  
  # average of communalities
  R2_aux <- R2[endo == 1]
  comu_aux <- n_comu <- NULL
  for (j in seq_len(lvs))
  {
    # non mono factorial blocks only
    if (sum(blocklist==j) > 1)
    {
      comu_aux = c(comu_aux, mean(comu[blocklist==j]))
      n_comu = c(n_comu, sum(blocklist==j))
    }
  }
  mean_communality = sum(comu_aux * n_comu) / sum(n_comu)
  gof = sqrt(mean_communality * mean(R2_aux))
  # output
  gof
}

get_rank <- function(X) 
{
  X_ranked = rep(0, length(X))
  uniq = unique(X, ties='min')
  uniq_ranked = rank(uniq)
  for (k in seq_len(length(uniq))) {
    X_ranked[which(X == uniq[k])] <- uniq_ranked[k]  
  }
  
  X_ranked
}

get_weights <- function(X, path_matrix, blocks, specs)
{
  lvs = nrow(path_matrix)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)-1) / nrow(X))   # std.dev factor correction
  blockinds = indexify(blocks)
  
  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  W = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w_old = rowSums(W)    
  iter = 1
  
  repeat 
  {
    # external estimation of LVs 'Y'
    Y = X %*% W
    Y = scale(Y) * sdv
    # matrix of inner weights 'e' 
    E <- switch(specs$scheme, 
                "centroid" = sign(cor(Y) * (path_matrix + t(path_matrix))),
                "factorial" = cor(Y) * (path_matrix + t(path_matrix)),
                "path" = get_path_scheme(path_matrix, Y))
    # internal estimation of LVs 'Z'
    Z = Y %*% E  
    #    Z = Z %*% diag(1/(apply(Z,2,sd)*sdv), lvs, lvs)  # scaling Z
    # computing outer weights 'w'
    for (j in seq_len(lvs))
    {
      if (specs$modes[j] == "A")
        W[blockinds==j,j] <- (1/nrow(X)) * Z[,j] %*% X[,blockinds==j] 
      if (specs$modes[j] == "B")
        W[blockinds==j,j] <- solve.qr(qr(X[,blockinds==j]), Z[,j])
    }
    w_new = rowSums(W)                
    w_dif = sum((abs(w_old) - abs(w_new))^2) 
    if (w_dif < specs$tol || iter == specs$maxiter) break
    w_old = w_new
    iter = iter + 1
  } # end repeat       
  
  # preparing results
  if (iter == specs$maxiter) {
    results = NULL
  } else {
    W = W %*% diag(1/(apply(X %*% W, 2, sd)*sdv), lvs, lvs)
    w_new = rowSums(W)                
    names(w_new) = colnames(X)
    dimnames(W) = list(colnames(X), rownames(path_matrix))    
    results = list(w = w_new, W = W, ODM = ODM, iter = iter)
  }
  # output
  results
}

is_numeric_vector <- function(x) {
  if (is.vector(x) & is.numeric(x)) TRUE else FALSE
}

is_numeric_dataframe <- function(x) {
  if (is.data.frame(x)) {
    numerics = unlist(lapply(x, is.numeric))
    if (sum(numerics) == dim(x)[2L]) TRUE else FALSE
  } else FALSE
}

is_string_dataframe <- function(x) {
  if (is.data.frame(x)) {
    characters = unlist(lapply(x, is.character))
    if (sum(characters) == dim(x)[2L]) TRUE else FALSE
  } else FALSE
}

is_string_vector <- function(x) {
  if (is.vector(x) & is.character(x)) TRUE else FALSE
}

is_logical_vector <- function(x) {
  if (is.vector(x) & is.logical(x)) TRUE else FALSE
}

is_not_vector <- function(x) {
  !is_vector(x)
}

get_weights_nonmetric <-
  function(X, path_matrix, blocks, specs)
  {
    lvs = nrow(path_matrix)
    mvs = ncol(X)
    num_obs = nrow(X)
    correction = sqrt(nrow(X) / (nrow(X)-1))
    blockinds = indexify(blocks)
    block_sizes = lengths(blocks)
    PLScomp = specs$plscomp
    start_end = from_to(block_sizes)
    
    # create dummy matrices for categorical manifest variables
    dummies = get_dummies(X, specs)
    
    ## transforming X in a list of blocks
    Xblocks = vector("list", lvs)
    start_end = from_to(blocks)
    from = start_end$from
    to = start_end$to
    for (q in seq_len(lvs)) {
      if (from[q] == to[q]) {
        Xblocks[[q]] = as.matrix(X[,from[q]:to[q]])
      } else {
        Xblocks[[q]] = X[,from[q]:to[q]]      
      }
    }
    
    # list for quantification of blocks' variables
    QQ = Xblocks
    # missing data flags
    missing_data = vapply(Xblocks, is_missing)
    # initialize list with availibility indicators (to handle NAs)
    X_avail = vector("list", lvs)
    # outer design matrix 'ODM' and matrix of outer weights 'W'
    ODM = list_to_dummy(blocks)
    
    # =======================================================================
    # initialization
    # =======================================================================
    # outer weights (normalized)
    w_ones = list_ones(block_sizes)
    w = lapply(w_ones, normalize)
    # LV scores
    Y = matrix(0, num_obs, lvs)
    for (q in seq_len(lvs)) {
      if (missing_data[q]) {
        # binary matrix (1=available data, 0=NA)
        X_avail[[q]] = 1 - is.na(Xblocks[[q]])
        for (i in seq_len(nrow(X))) {
          aux_numerator = sum(QQ[[q]][i,]*w[[q]], na.rm = TRUE)
          aux_denom = sum(w[[q]][which(is.na(QQ[[q]][i,]*w[[q]]) == FALSE)]^2)
          Y[i,q] <- aux_numerator / aux_denom
        }
      } else {
        Y[,q] = QQ[[q]] %*% w[[q]]        
      }
    }
    
    # =======================================================================
    # iterative cycle
    # =======================================================================
    # matrix of inner weights
    E = matrix(0, lvs, lvs)
    link = t(path_matrix) + path_matrix
    z_temp = matrix(0, num_obs, 1)
    iter = 0
    repeat 
    {
      iter = iter + 1
      #    y_old = as.numeric(Y)
      Y_old = Y
      
      # =============================================================
      # updating inner weights
      # =============================================================
      E <- switch(specs$scheme, 
                  "centroid" = sign(cor(Y) * link),
                  "factorial" = cov(Y) * link,
                  "path" = get_path_scheme(path_matrix, Y))
      # internal estimation of LVs 'Z'
      Z = Y %*% E
      
      # for each block
      for (q in seq_len(lvs))
      {
        # standardize inner estimates if PLScore mode
        # if (specs$modes[q] != "PLSCOW") {
        #### Giorgio's suggestion: do not standardize inner estimates:
        #if (specs$modes[q] != "PLSCOW" & specs$modes[q] != "NEWA") {
        #  Z[,q] <- scale(Z[,q]) * correction
        #}
        # =============================================================
        # Quantification of the MVs in block ["QQ"]
        # =============================================================
        # for each MV in block 'q'
        if (specs$modes[q] == "B" && block_sizes[q] > 1) {
          Beta <- summary(lm(Z[,q]~QQ[[q]]))$coef[-1,1]
          X.star <- matrix(,num_obs,block_sizes[q])
          for (p in 1L:block_sizes[q]) {
            X.star[,p] <- (1/Beta[p])*(Z[,q] - (QQ[[q]][,-p,drop=FALSE]%*%Beta[-p]))
            if (specs$scaling[[q]][p] == "nom") {
              # extract corresponding dummy matrix
              #          aux_dummy = dummies[[blocks[[q]][p]]]
              which_dummy = (start_end$from[q]:start_end$to[q])[p]
              aux_dummy = dummies[[which_dummy]]
              # apply scaling
              QQ[[q]][,p] = get_nom_scale(X.star[,p], Xblocks[[q]][,p], aux_dummy)
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }
            if (specs$scaling[[q]][p] == "ord") {
              # extract corresponding dummy matrix
              #          aux_dummy = dummies[[blocks[[q]][p]]]
              which_dummy = (start_end$from[q]:start_end$to[q])[p]
              aux_dummy = dummies[[which_dummy]]
              # apply scaling
              QQ[[q]][,p] = get_ord_scale(X.star[,p], Xblocks[[q]][,p], aux_dummy)
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }                   
            if (specs$scaling[[q]][p] == "num") {
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }
          }
        }
        else {
          for (p in 1L:block_sizes[q]) {
            if (specs$scaling[[q]][p] == "nom") {
              # extract corresponding dummy matrix
              #          aux_dummy = dummies[[blocks[[q]][p]]]
              which_dummy = (start_end$from[q]:start_end$to[q])[p]
              aux_dummy = dummies[[which_dummy]]
              # apply scaling
              QQ[[q]][,p] = get_nom_scale(Z[,q], Xblocks[[q]][,p], aux_dummy)
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }
            if (specs$scaling[[q]][p] == "ord") {
              # extract corresponding dummy matrix
              #          aux_dummy = dummies[[blocks[[q]][p]]]
              which_dummy = (start_end$from[q]:start_end$to[q])[p]
              aux_dummy = dummies[[which_dummy]]
              # apply scaling
              QQ[[q]][,p] = get_ord_scale(Z[,q], Xblocks[[q]][,p], aux_dummy)
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }                   
            if (specs$scaling[[q]][p] == "num") {
              QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
            }
            ### DO WE REALLY NEED THIS LINE:
            #if (specs$scaling[[q]][p] == "raw") {
            #  QQ[[q]][,p] = QQ[[q]][,p]
            #}
          }
        }
        
        # =============================================================
        # updating outer weights "w" and outer estimates "Y"
        # =============================================================
        
        # Mode A (="PLScore1comp") ====================================
        if (specs$modes[q] == "A") {
          if (missing_data[q]) {
            # compute w[[q]][l] as the regr. coeff. of QQ[[q]][,l] on Z[,q] 
            # considering only the lines where QQ[[q]][i,l] exist
            w[[q]] = colSums(QQ[[q]]*Z[,q], na.rm = TRUE)
            w[[q]] = w[[q]] / colSums((X_avail[[q]]*Z[,q])^2)
            # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
            # considering only the columns where QQ[[q]][i,l] exist
            Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm = TRUE)
            Y[,q] = Y[,q] / colSums((t(X_avail[[q]])*w[[q]])^2)
            # normalize Y[,q] to unitary variance
            Y[,q] = scale(Y[,q]) * correction
          } 
          else {# complete data in block q
            w[[q]] = (t(QQ[[q]]) %*% Z[,q]) / sum(Z[,q]^2)
            Y[,q] = QQ[[q]] %*% w[[q]]
            Y[,q] = scale(Y[,q]) * correction
          }
        }
        
        #  Mode New A (= "PLScow1comp") ================================
        if (specs$modes[q] == "NEWA") {
          if (missing_data[q]) {
            # compute w[[q]][l] as the regr. coeff. of QQ[[q]][,l] on Z[,q]
            # considering just the lines where QQ[[q]][i,l] exist
            w[[q]] = colSums(QQ[[q]]*Z[,q], na.rm = TRUE)
            w[[q]] = w[[q]]/colSums((X_avail[[q]]*Z[,q])^2)
            # normalize w[[q]] to unitary norm
            w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
            # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
            # considering only the columns where QQ[[q]][i,l] exist
            Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
            Y[,q] = Y[,q] / colSums((t(X_avail[[q]])*w[[q]])^2)
          } 
          else {# complete data in block q
            w[[q]] = (t(QQ[[q]]) %*% Z[,q]) / sum(Z[,q]^2)
            w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
            Y[,q] = QQ[[q]] %*% w[[q]]
          }
        }
        
        # Mode B (NAs were not allowed.. so far. Now we can use PLSR) ====
        if (specs$modes[q] == "B") {
          if (missing_data[q]) {# use full component PLS-R 
            w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = block_sizes[q])$B
            # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
            # considering only the columns where QQ[[q]][i,l] exist
            Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
            Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
            # normalize Y[,q] to unitary variance
            Y[,q] = scale(Y[,q]) * correction		
          }
          else {# complete data in block q
            w[[q]] = solve.qr(qr(QQ[[q]]), Z[,q])
            #w[[q]] = solve(t(QQ[[q]]) %*% QQ[[q]]) %*% t(QQ[[q]]) %*% Z[,q]
            Y[,q] = QQ[[q]] %*% w[[q]]
            Y[,q] = scale(Y[,q]) * correction
          }	
        }
        
        # Mode PLScore ===================================================
        if (specs$modes[q] == "PLSCORE") {
          if (missing_data[q]) {
            w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
            # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
            # considering only the columns where QQ[[q]][i,l] exist
            Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
            Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
            # normalize Y[,q] to unitary variance
            Y[,q] = scale(Y[,q]) * correction		
          }
          else {# complete data in block q
            w[[q]] = get_PLSR(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
            Y[,q] = QQ[[q]] %*% w[[q]]
            Y[,q] = scale(Y[,q]) * correction
          }	
        }
        
        # Mode PLScow =====================================================
        if (specs$modes[q] == "PLSCOW") {
          if (missing_data[q]) {
            w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
            # normalize w[[q]] to unitary norm
            w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
            # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
            # considering only the columns where QQ[[q]][i,l] exist
            Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
            Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
          }
          else {# complete data in block q
            w[[q]] = get_PLSR(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
            w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
            Y[,q] = QQ[[q]] %*% w[[q]]
          }	
        }
        
      }
      # check convergence
      convergence <- sum((abs(Y_old) - abs(Y))^2)
      # Y_old: keep it as a matrix
      #    convergence <- sum((abs(y_old) - abs(as.numeric(Y)))^2)
      if (convergence < specs$tol | iter > specs$maxiter) 
        break
    } # end repeat
    
    # preparing results
    if (iter == specs$maxiter) {
      results = NULL
    } else {
      W = list_to_matrix(lapply(w, as.numeric))
      # open new lines
      QQ = do.call("cbind", QQ)
      W = W %*% diag(1/(apply(QQ %*% W, 2, sd, na.rm=TRUE)/correction), lvs, lvs)
      # end new lines
      w = rowSums(W)
      dimnames(W) = list(colnames(X), rownames(path_matrix))
      dimnames(Y) = list(rownames(X), rownames(path_matrix))
      results = list(w = w, W = W, Y = Y, QQ = QQ, ODM = ODM, iter = iter)
    }
    # output
    results  
  }


get_boot_stats <- function(MATRIX, original) {
  data.frame(Original = original,
             Mean.Boot = apply(MATRIX, 2, mean), 
             Std.Error = apply(MATRIX, 2, sd), 
             perc.025 = apply(MATRIX, 2, function(x) quantile(x, 0.025)),
             perc.975 = apply(MATRIX, 2, function(x) quantile(x, 0.975)))
}

from_to <- function(x, ...) {
  UseMethod("from_to", x)  
}


#' @S3method from_to default
from_to.default <- function(x, ...)
{
  if (!is_numeric_vector(x) || !list_of_vectors(x))
    stop("\n'from_to()' requires a numeric vector or a list of vectors")
}


#' @S3method from_to numeric
from_to.numeric <- function(x, ...)
{
  if (!is_numeric_vector(x))
    stop("\n'from_to()' requires a numeric vector")
  
  to = cumsum(x)
  from = to - x + 1
  list(from=from, to=to)
}


#' @S3method from_to list
from_to.list <- function(x, ...)
{
  if (!list_of_vectors(x))
    stop("\n'from_to()' requires a list of vectors")
  
  aux = unlist(lapply(x, length))
  to = cumsum(aux)
  from = to - aux + 1
  list(from=from, to=to)
}

get_path_scheme <- function(path_matrix, LV)
{
  # matrix for inner weights
  E = path_matrix
  
  for (k in seq_len(ncol(path_matrix))) 
  {
    # followers
    follow <- path_matrix[k,] == 1
    if (sum(follow) > 0)
      E[follow,k] <- lm(LV[,k] ~ LV[,follow] - 1)$coef
    # predecesors
    predec <- path_matrix[,k] == 1
    if (sum(predec) > 0)
      E[predec,k] <- cor(LV[,k], LV[,predec])
  } 
  
  # output
  E
}

get_dummies <- function(MV, specs)
{  
  # get metric
  metric = get_metric(specs$scaling)
  
  if (metric) {
    dummies = NULL
  } else {
    scaling = unlist(specs$scaling)
    nominal_ordinal = which(scaling %in% c("ord", "nom"))
    dummies = vector(mode="list", length(scaling))
    # only categorical mvs have an associated dummy matrix
    for (j in seq_along(nominal_ordinal)) {
      dummies[[nominal_ordinal[j]]] = get_dummy(MV[,nominal_ordinal[j]])
    }    
  }
  
  # output
  dummies
}

get_nom_scale <- function(y, x, Xdummy) 
{
  n = length(x)
  p = max(x, na.rm = TRUE)
  
  # vector of the means of y conditioned to the levels of x
  quant = tapply(y, x, mean, na.rm=TRUE)
  x_quant = Xdummy %*% quant
  x_quant
  
  # ===========  just in the case you need of them ============
  # ===========  eta2 = correlation rato  ============
  #eta2 <- var(x_quant)/var(y)
  #eta2<-var(x_quant,na.rm = T)/var(y, na.rm = T)	
  #list(Xdummy=Xdummy,x_quant=x_quant, eta2=eta2, quant=quant)
}

get_num_scale <- function(X) {
  X = as.matrix(X)
  X_scaled = matrix(0, nrow(X), ncol(X))
  for (j in seq_len(ncol(X) )) {
    correction <- (sqrt(length(na.omit(X[,j]))/(length(na.omit(X[,j]))-1)))
    X_scaled[,j] <- scale(X[,j]) * correction
  }
  #rownames(X_scaled) = rownames(X) 
  #colnames(X_scaled) = colnames(X)
  X_scaled
}

get_ord_scale <- function(y, x, Xdummy) 
{
  n <- length(x)
  p <- max(x, na.rm = TRUE)
  #	# ===========  build the (p x p) matrix of zeros ====	
  #	Xdummy<-matrix(0,n,p)
  #	# ===========  put the ones ============
  #	for (k in 1:p) {
  #		Xdummy[x == k,k] = 1
  #  	}
  #  	# ===========  if there are NA, add them ============
  #	if (any(is.na(x))) {
  #		Xdummy[which(rowSums(Xdummy) == 0),] <- NA
  #	}
  # =====  building an initial vector of scaling parameters ======
  quant <- (tapply(y, x, mean, na.rm=TRUE))
  
  # =====  searching for monotonically increasing quantifications ======
  quant_incr <- quant
  Xdummy_incr <- Xdummy
  repeat {
    ncol_Xdummy_Old <- ncol(Xdummy_incr)
    # ===== if the monotony is not respected, merge the columns =====	
    for (k in seq_len(ncol(Xdummy_incr)-1)) {
      if (quant_incr[k] > quant_incr[k+1]) {
        Xdummy_incr[,k+1] <- Xdummy_incr[,k] + Xdummy_incr[,k+1]
        Xdummy_incr <- as.matrix(Xdummy_incr[,-k])
        quant_incr <- c()
        for (k in seq_len(ncol(Xdummy_incr))) {
          quant_incr[k] <- sum((Xdummy_incr[,k])*y,na.rm=TRUE)
          quant_incr[k] <- quant_incr[k]/sum(Xdummy_incr[,k], na.rm=TRUE)
        }
        break
      }
    }
    if (ncol(Xdummy_incr) == 1 || (ncol(Xdummy_incr) == ncol_Xdummy_Old)) {break}
  }
  x_quant_incr <- Xdummy_incr %*% quant_incr
  
  var_incr <- var(x_quant_incr, na.rm = TRUE)
  
  # =====  searching for monotonically decreasing quantifications ======
  quant_decr <- quant
  Xdummy_decr <- Xdummy
  repeat {
    ncol_Xdummy_Old<-ncol(Xdummy_decr)
    # ===== if the monotony is not respected, merge the columns =====		
    for (k in seq_len(ncol(Xdummy_decr)-1)) {
      if (quant_decr[k] < quant_decr[k+1]) {
        Xdummy_decr[,k+1] <- Xdummy_decr[,k] + Xdummy_decr[,k+1]
        Xdummy_decr <- as.matrix(Xdummy_decr[,-k])
        quant_decr <- c()
        for (k in seq_len(ncol(Xdummy_decr))) {
          quant_decr[k] <- sum((Xdummy_decr[,k])*y,na.rm=TRUE)
          quant_decr[k] <- quant_decr[k]/sum(Xdummy_decr[,k], na.rm=TRUE)
        }
        break
      }
    }
    if (ncol(Xdummy_decr) == 1 || (ncol(Xdummy_decr) == ncol_Xdummy_Old)) {break}
  }
  x_quant_decr <- Xdummy_decr %*% quant_decr
  var_decr <- var(x_quant_decr, na.rm = TRUE)
  
  # =====  choosing between increasing and decreasing ======
  if (var_incr < var_decr) {
    Xdummy <- Xdummy_decr
    quant <- -(quant_decr)	
    x_quant <- -(x_quant_decr)	
  }
  else {
    Xdummy <- Xdummy_incr
    quant <- quant_incr	
    x_quant <- x_quant_incr	
  }
  
  x_quant
  # ===========  just in the case you need of them ============
  # ===========  eta2 = correlation rato  ============
  #eta2 <- var(x_quant)/var(y)
  #eta2<-var(x_quant, na.rm = T)/var(y, na.rm = T)	
  #list(Xdummy=Xdummy,x_quant=x_quant, eta2=eta2, quant=quant)
}

get_PLSR <- function(Y, X, ncomp) {
  Y = as.matrix(Y)
  X = as.matrix(X)
  colnamesX = colnames(X)
  n = nrow(Y)
  p = ncol(X)
  A = matrix(NA, p, ncomp)
  rownames(A) <- colnamesX
  colnames(A) <- c(seq_len(ncomp))
  T <- matrix(NA, n, ncomp)
  rownames(T) <- rownames(X)
  colnames(T) <- c(seq_len(ncomp))
  C <- matrix(NA, 1, ncomp)
  rownames(C) <- c("z")
  colnames(C) <- c(seq_len(ncomp))
  P <- matrix(NA, p, ncomp)
  rownames(P) <- colnamesX
  colnames(P) <- c(seq_len(ncomp))
  Wstar <- matrix(NA, p, ncomp)
  B <- matrix(,p,1)
  rownames(B) <- colnamesX
  colnames(B) <- c("z")
  R2 <- matrix(NA, 2, ncomp)
  varX <- sum(apply(X,2,var))
  varY <- var(Y)
  for (k in seq_len(ncomp)) {
    if (k == 1) {
      A[,k] <- t(X) %*% Y
      A[,k] <- A[,k] / sqrt(sum(A[,k]^2))
      T[,k] <- X %*% A[,k]
      C[1,k] <- t(Y) %*% T[,k]/(sum(T[,k]^2))
      P[,k] <- t(X) %*% T[,k]/(sum(T[,k]^2))
      Xres <- X - T[,k] %*% t(P[,k])		
      if (ncomp == 1) {
        B <- A[,k] * C[1,k]
      }
      Wstar <- as.matrix(A[,k])
      R2[1,k] <- 1 - (sum(apply(Xres,2,var)) / varX)
      R2[2,k] <- 1 - (var(Y-(T[,k]%*%t(C[,k]))) / varY)
    }
    else {
      A[,k] <- t(Xres) %*% Y
      A[,k] <- A[,k] / sqrt(sum(A[,k]^2))
      T[,k] <- Xres %*% A[,k]
      C[1,k] <- t(Y) %*% T[,k]/(sum(T[,k]^2))
      P[,k] <- t(Xres) %*% T[,k]/(sum(T[,k]^2))
      Xres <- Xres - T[,k] %*% t(P[,k])
      R2[1,k] <- 1 - (sum(apply(Xres,2,var)) / varX)
      R2[2,k] <- 1 - (var(Y-(T[,seq_len(k)] %*% as.matrix(C[,seq_len(k)]))) / varY)
    }
  }
  Wstar <- A %*% solve(t(P) %*% A)	
  B <- Wstar %*% t(C)
  Pcorr <- as.matrix(cor(X,T))
  Ccorr <- as.matrix(cor(Y,T))
  rownames(R2) <- c("X","z")
  colnames(R2) <- c(seq_len(ncomp))
  rownames(Pcorr) <- colnamesX
  colnames(Pcorr) <- c(seq_len(ncomp))
  rownames(Ccorr) <- c("z")
  colnames(Ccorr) <- c(seq_len(ncomp))
  # output
  list(T = T, 
       C = C, 
       P = P, 
       A = A, 
       B = B, 
       R2cum = R2, 
       Pcorr = Pcorr, 
       Ccorr = Ccorr, 
       Wstar = Wstar)
}

get_PLSR_NA <- function(Y, X, ncomp) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  X_avail <- 1-is.na(X)
  colnamesX <- colnames(X)
  n <- nrow(Y)
  p <- ncol(X)
  A <- matrix(,p,ncomp)
  rownames(A) <- colnamesX
  colnames(A) <- c(seq_len(ncomp))
  T <- matrix(,n,ncomp)
  rownames(T) <- rownames(X)
  colnames(T) <- c(seq_len(ncomp))
  C <- matrix(,1,ncomp)
  rownames(C) <- c("z")
  colnames(C) <- c(seq_len(ncomp))
  P <- matrix(,p,ncomp)
  rownames(P) <- colnamesX
  colnames(P) <- c(seq_len(ncomp))
  Wstar <- matrix(,p,ncomp)
  B <- matrix(,p,1)
  rownames(B) <- colnamesX
  colnames(B) <- c("z")
  R2 <- matrix(,2,ncomp)
  varX <- sum(apply(X,2,function(x){var(x,na.rm=TRUE)}))
  varY <- var(Y)
  for (k in seq_len(ncomp)) {
    if (k == 1) {
      A[,k] <- colSums(X*Y[,1], na.rm = TRUE)
      A[,k] <- A[,k]/colSums((X_avail*Y[,1])^2) 
      A[,k] <- A[,k]/sqrt(sum(A[,k]^2))
      T[,k] <- colSums(t(X)*A[,k], na.rm=TRUE) 
      T[,k] <- T[,k]/colSums((t(X_avail)*A[,k])^2)
      C[1,k] <- sum(Y[,1]*T[,k], na.rm=TRUE)/(sum(T[,k]^2, na.rm=TRUE))
      P[,k] <- colSums(X*T[,k],na.rm=TRUE)
      P[,k] <- P[,k]/colSums((X_avail*T[,k])^2) 
      Xres <- X - T[,k]%*%t(P[,k])
      Yres <- Y - T[,k]*C[1,k]		
      if (ncomp == 1) {
        B <- A[,k]*C[1,k]
      }
      Wstar <- as.matrix(A[,k])
      R2[1,k] <-1 - (sum(apply(Xres,2,function(x){var(x,na.rm=TRUE)})) / varX)
      R2[2,k] <- 1 - (var(Yres)/varY)
    }
    else {
      A[,k] <-colSums(Xres*Yres[,1],na.rm=TRUE)
      A[,k] <- A[,k]/colSums((X_avail*Yres[,1])^2) 
      A[,k] <- A[,k]/sqrt(sum(A[,k]^2))
      T[,k] <- colSums(t(Xres)*A[,k], na.rm=TRUE) 
      T[,k] <- T[,k]/colSums((t(X_avail)*A[,k])^2)
      C[1,k] <- sum(Yres[,1]*T[,k], na.rm=TRUE)/(sum(T[,k]^2, na.rm=TRUE))
      P[,k] <- colSums(Xres*T[,k],na.rm=TRUE)  
      P[,k] <- P[,k]/colSums((X_avail*T[,k])^2)
      Xres <- Xres - T[,k]%*%t(P[,k])
      Yres <- Yres - T[,k]*C[1,k]
      R2[1,k] <-1 - (sum(apply(Xres,2,function(x){var(x,na.rm=TRUE)})) / varX)
      R2[2,k] <- 1 - (var(Yres)/varY)
    }
  }
  uppertri_PA <- ((t(P)%*%A)*upper.tri(diag(ncomp)))+diag(ncomp)
  Wstar <- A%*%solve(uppertri_PA)	
  B <- Wstar%*%t(C)
  Pcorr <- as.matrix(cor(X,T,use="pairwise.complete.obs"))
  Ccorr <- as.matrix(cor(Y,T))
  rownames(Wstar) <- colnamesX
  colnames(Wstar) <- c(seq_len(ncomp))
  rownames(R2) <- c("X","z")
  colnames(R2) <- c(seq_len(ncomp))
  rownames(Pcorr) <- colnamesX
  colnames(Pcorr) <- c(seq_len(ncomp))
  rownames(Ccorr) <- c("z")
  colnames(Ccorr) <- c(seq_len(ncomp))
  # output
  list(T = T, 
       C = C, 
       P = P, 
       A = A, 
       B = B, 
       R2cum = R2, 
       Pcorr = Pcorr, 
       Ccorr = Ccorr, 
       Wstar = Wstar, 
       Xres = Xres, 
       Yres = Yres)
}

get_dummy <- function(x) 
{
  n = length(x)
  # p = max(x, na.rm = TRUE)
  # since 'x' could include zero, it's better to use the following:
  p = length(unique(x[!is.na(x)]))
  
  # build the (p x p) dummy matrix 
  Xdummy = matrix(0, n, p)
  for (k in seq_len(p)) {
    Xdummy[x == k,k] = 1
  }
  # if there are NAs, add them
  if (any(is.na(x))) {
    Xdummy[which(rowSums(Xdummy) == 0),] <- NA
  }
  # output
  Xdummy
}

get_boots <-
  function(DM, path_matrix, blocks, specs, br)
  {
    # =======================================================
    # inputs setting
    # =======================================================  
    lvs = nrow(path_matrix)
    lvs.names = rownames(path_matrix)
    mvs = ncol(DM)
    mvs.names = colnames(DM)
    blocklist = indexify(blocks)
    endo = sign(rowSums(path_matrix))
    bootnum = br
    # apply corresponding treatment (centering, reducing, ranking)
    X = get_treated_data(DM, specs)
    
    # =======================================================
    # computation of the original plspm model
    # =======================================================  
    metric = get_metric(specs$scaling)
    if (metric) {
      # object 'weights' contains outer w's, W, ODM, iter
      out.ws = get_weights(X, path_matrix, blocks, specs)
      ok_weights = test_null_weights(out.ws, specs)
      wgs.orig = out.ws$w
      Y.lvs = get_scores(X, out.ws$W)
    } else {
      # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
      out.ws = get_weights_nonmetric(X, path_matrix, blocks, specs)
      ok_weights = test_null_weights(out.ws, specs)
      wgs.orig = out.ws$w
      Y.lvs = out.ws$Y
      X = out.ws$QQ  # quantified MVs
    }
    
    pathmod <- get_paths(path_matrix, Y.lvs)
    Path <- pathmod[[2]]
    path.orig <- as.vector(Path[path_matrix==1])
    r2.orig <- pathmod[[3]][endo==1]
    Path.efs <- get_effects(Path)
    xloads = cor(X, Y.lvs)
    loads.orig = rowSums(xloads * out.ws$ODM)
    
    # =======================================================
    # Bootstrap Validation
    # =======================================================  
    path.labs <- NULL
    efs.labs <- NULL
    for (j in seq_len(lvs))
      for (i in j:lvs)
        if (path_matrix[i,j]==1) 
          path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i]))
    WEIGS <- matrix(NA, bootnum, mvs)
    LOADS <- matrix(NA, bootnum, mvs)
    PATHS <- matrix(NA, bootnum, sum(path_matrix))
    TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
    RSQRS <- matrix(NA, bootnum, sum(endo))
    i <- 1
    while (i <= bootnum)
    {
      boot.obs <- sample.int(nrow(X), size=nrow(X), replace=TRUE)
      DM.boot <- DM[boot.obs,]
      # apply corresponding treatment (centering, reducing, ranking)
      X.boot = get_treated_data(DM.boot, specs)        
      # calculating boot model parameters 
      if (metric) {
        # object 'weights' contains outer w's, W, ODM, iter
        w.boot = get_weights(X.boot, path_matrix, blocks, specs)
        if (is.null(w.boot)) {
          i <- i - 1
          next
        }
        Y.boot = get_scores(X.boot, w.boot$W)
      } else {
        # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
        w.boot = get_weights_nonmetric(X.boot, path_matrix, blocks, specs)
        if (is.null(w.boot)) {
          i <- i - 1
          next
        }
        Y.boot = w.boot$Y
        X.boot = w.boot$QQ  # quantified MVs
        # X.boot = do.call(cbind, w.boot$QQ)  # quantified MVs
      }
      WEIGS[i,] <- w.boot$w
      pathmod <- get_paths(path_matrix, Y.boot)
      P.boot <- pathmod[[2]]
      Toef.boot <- get_effects(P.boot)
      PATHS[i,] <- as.vector(P.boot[path_matrix==1])
      TOEFS[i,] <- Toef.boot[,4]
      RSQRS[i,] <- pathmod[[3]][endo==1]
      xloads = cor(X.boot, Y.boot)
      LOADS[i,] = rowSums(xloads * w.boot$ODM)
      i <- i + 1
    }
    
    # =======================================================
    # Bootstrap results
    # ======================================================= 
    # Outer weights
    WB = get_boot_stats(WEIGS, wgs.orig)
    #rownames(WB) = mvs.names
    rownames(WB) <- paste(rep(lvs.names, vapply(blocks, length)), 
                          mvs.names, sep='-')
    # Loadings
    LB = get_boot_stats(LOADS, loads.orig)
    #rownames(LB) = mvs.names
    rownames(LB) <- paste(rep(lvs.names, vapply(blocks, length)), 
                          mvs.names, sep='-')
    # Path coefficients
    colnames(PATHS) = path.labs
    PB = get_boot_stats(PATHS, path.orig)
    # R-squared
    colnames(RSQRS) = lvs.names[endo == 1]
    RB = get_boot_stats(RSQRS, r2.orig)
    # Total effects
    colnames(TOEFS) = Path.efs[, 1]
    TE = get_boot_stats(TOEFS, Path.efs[,4]) 
    
    # Bootstrap Results
    list(weights = WB, 
         loadings = LB, 
         paths = PB, 
         rsq = RB, 
         total.efs = TE)
  }

list_ones <- function(x)
{
  if (!is_numeric_vector(x))
    stop("\n'list_ones()' requires a numeric vector")
  # output
  lapply(x, function(u) rep(1, u))
}

list_to_matrix <- function(alist)
{
  if (!list_of_numeric_vectors(alist))
    stop("\n'list_to_matrix()' requires a list of numeric vectors")
  
  aux = lengths(alist)
  to = cumsum(aux)
  from = to - aux + 1
  linked_matrix = matrix(0, sum(aux), length(aux))
  for (j in seq_along(aux))
    linked_matrix[from[j]:to[j], j] = alist[[j]]
  linked_matrix
}

is_vector <- function(x) {
  if (is.vector(x)) TRUE else FALSE
}
