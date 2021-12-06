# Helper Functions --------------------------------------------------------

# Simple function for passing list arguments to mapply
papply <- function(.l, .f, ...) {
  args <- c(.l, list(FUN = .f, MoreArgs = list(...), SIMPLIFY = FALSE))
  return(do.call(mapply, args = args))
}

get_HLSM_type <- function(object_list) {
  calls <- lapply(object_list, getCall)
  funcs <- sapply(lapply(calls, `[[`, 1), as.character)
  type <- unique(gsub('HLSM(.*)EF', '\\1', funcs))
  if (length(type) > 1) {
    stop("HLSM list must be all of the same type.")
  } else if (!(type %in% c('fixed', 'random', 'LSM'))) {
    stop("Unknown HLSM type found in object.")
  }
  if(type=="LSM"){
  	test=deparse(calls[[1]]) #creates string of the call
  	est.int=grep("estimate.intercept = TRUE", test)
  	if(length(est.int)>0){type="LSM.estInt"}else{type="LSM.fixedInt"}
  }
  return(type)
}


# Draws Extraction & Conversion --------------------------------------------

extract_param <- function(chain, type, burnin = 0, thin = 1) {
  # in random, niter X nnet matrix
  
  if(type=="LSM.fixedInt"){
  	 beta_draws <- getBeta(chain, burnin = burnin, thin = thin)
  
  # Creating shape and dimnames to pass to array creation, with goal of
  # binding intercept to beta array along the "variable" axis.
  beta_shape <- dim(beta_draws)
  beta_dnames <- list(
    iterations = seq_len(beta_shape[1]),
    variables = paste0('X', seq_len(beta_shape[2]))
  )}else{
  inter_draws <- getIntercept(chain, burnin = burnin, thin = thin)
  # in random, niter X nvar X nnet matrix
  beta_draws <- getBeta(chain, burnin = burnin, thin = thin)
  
  # Creating shape and dimnames to pass to array creation, with goal of
  # binding intercept to beta array along the "variable" axis.
  beta_shape <- dim(beta_draws)
  beta_dnames <- list(
    iterations = seq_len(beta_shape[1]),
    variables = paste0('X', seq_len(beta_shape[2]))
  )
  inter_shape <- beta_shape
  inter_shape[2] <- 1
  inter_dnames <- list(
    iterations = seq_len(inter_shape[1]),
    variables = 'Intercept'
  )}
  
  # If model is not random effects, it is fixed across network
  if (type == "random") {
    beta_dnames$network <- paste0('Net', seq_len(beta_shape[3]))
    inter_dnames$network <- paste0('Net', seq_len(inter_shape[3]))
  } else if (type != "fixed" && type != "LSM.estInt" && type !="LSM.fixedInt") {
    stop("Type must be either 'LSM', 'fixed', or 'random'")
  }
  
  # Apply dnames and new variable dimension to intercept array, then bind along
  # the variable dimension.
  beta_array <- array(beta_draws, dim = beta_shape, dimnames = beta_dnames)
  if(type=="LSM.fixedInt"){
  	param_array=beta_array
  }else{
  inter_array <- array(inter_draws, dim = inter_shape, dimnames = inter_dnames)
  param_array <- abind(inter_array, beta_array, along = 2)
  }
  return(param_array)
}

param_to_mcmc <- function(param) {
  param_df <- as.data.frame(param)
  param_mcmc <- as.mcmc(param_df)
  return(param_mcmc)
}


# PSRF Functions ----------------------------------------------------------

psrf_param <- function(param_mcmc_list, warn_1chain = TRUE) {
  result <- NULL
  if (nchain(param_mcmc_list) > 1) {
    result <- gelman.diag(param_mcmc_list, autoburnin = FALSE)
  } else if (warn_1chain) {
    warning("You have only provided one chain. PSRF is not available and will ",
            "not be included in the output.")
  }
  return(result)
}

psrf_summary <- function(result) {
  psrf_table <- result$psrf
  
  mask <- psrf_table[, 'Point est.'] > 1.05
  if (sum(mask) == 0) {
    mask <- which.max(psrf_table[, 'Upper C.I.'])
  }
  # if only 1 row, doesn't reduce to vector
  bad_psrf_table <- psrf_table[mask, , drop = FALSE]
  
  max_lim_psrf <- max(psrf_table[, 'Upper C.I.'])
  multi_psrf <- result$mpsrf
  
  cat("Potential Scale Reduction Factor:\n")
  cat("Gelman-Rubin between-chain convergence diagnostic.\n")
  cat("Upper C.I. near 1 indicates convergence.\n")
  cat("---\n")
  cat("Variable(s) with worst convergence\n")
  print(bad_psrf_table)
  cat("---\n")
  cat("Maximum Upper C.I.: ", round(max_lim_psrf, 4), '\n')
  cat("Multivariate PSRF Point Estimate: ", multi_psrf, '\n')
}

# Raftery Functions -------------------------------------------------------

raftery_param <- function(param_mcmc) {
  result <- NULL
  if (dim(param_mcmc)[1] <= 3746) {
    warning(
      "The chain length is less than the raftery diagnostic minimum length of ",
      "3746.\n",
      "If would like the raftery diagnostic information, ensure the chain ",
      "length of > 3746 iterations."
    )
  } else {
    result <- as.data.frame(raftery.diag(param_mcmc)$resmatrix)
    result$Nmin <- NULL
    colnames(result) <- c("burnin", "niters", "thinning")
  }
  return(result)
}

raftery_summary <- function(results) {
  longest_stats <- lapply(results, apply, 2, max)
  chain_stats <- do.call(rbind, longest_stats)
  rownames(chain_stats) <- paste("Chain", seq_along(results))
  cat("Raftery Diagnostics:\n")
  cat("Suggested Chain Specifications\n")
  print(chain_stats[, c(2, 1, 3)])
}


# Plotting Functions ------------------------------------------------------

plot_shape <- function(param_mcmc, type) {
  if (type == "fixed" | type == "LSM.fixedInt" | type == "LSM.estInt") {
    n_vars <- length(varnames(param_mcmc))
    nrows <- min(4, n_vars)
    ncols <- 1
    plotdex <- seq_len(n_vars)
  } else if (type == "random") {
    # Determine the number of distinct variables nad networks in the mcmc object
    vars <- varnames(param_mcmc)
    var_splits <- strsplit(vars, '.Net')
    uvars <- unique(lapply(var_splits, `[`, 1))
    unets <- unique(lapply(var_splits, `[`, 2))
    n_uvars <- length(uvars)
    n_unets <- length(unets)
    
    # set the number of columns and rows for the plot
    nrows <- min(4, n_unets)
    ncols <- min(3, n_uvars)
    
    # create a matrix of indices to control the plot order. this will plot each
    # variable in its own row. the networks will be spaced throughout all of the
    # networks
    big_dex_mat <- t(matrix(seq_len(n_uvars * n_unets), ncol = n_unets))
    netdex <- floor(seq(1, n_unets, length.out = ncols))
    dex_mat <- big_dex_mat[netdex,]
    plotdex <- as.vector(dex_mat)
  }else{
    stop("Type must be either 'fixed' or 'random'.")
  }
  
  return(list(nrows = nrows, ncols = ncols, plotdex = plotdex))
}

param_get_acf <- function(param_mcmc) {
  param_ts <- apply(param_mcmc, 2, as.ts)
  results <- apply(param_ts, 2, acf, plot = FALSE)
  lags <- lapply(results, `[[`, 'lag')
  acfs <- lapply(results, `[[`, 'acf')
  return(list(lag = lags, acf = acfs))
}

autocorr_param <- function(param_mcmc_list, col = 1:6, lty = 1) {
  vars <- varnames(param_mcmc_list)
  acf_results <- lapply(param_mcmc_list, param_get_acf)
  acf_results_t <- papply(acf_results, list)
  lags_bind <- papply(acf_results_t$lag, cbind)
  acf_bind <- papply(acf_results_t$acf, cbind)
  for (i in seq_len(nvar(param_mcmc_list))) {
    main <- paste("ACF Plot of", vars[i])
    nchains <- nchain(param_mcmc_list)
    plot_lags <- jitter(lags_bind[[i]], ifelse(nchains - 1, 2, 0))
    matplot(plot_lags, acf_bind[[i]], type = 'h', col = col, lty = lty,
            main = main)
  }
}


# Main Function -----------------------------------------------------------

HLSMdiag <- function(object, burnin = 0,
                     diags = c('psrf', 'raftery', 'traceplot', 'autocorr'),
                     col = 1:6, lty = 1) {
  if (is(object, 'HLSM')) {
    object_list <- list(object)
  } else if (is(object, 'list')) {
    object_list <- object
  } else {
    stop("object must be single HLSM chain or list of HLSM chains")
  }
  
  warn_1chain <- !missing(diags)
  # the default behavior is to return all information
  diags <- match.arg(diags, several.ok = TRUE)
  
  type <- get_HLSM_type(object_list)
  
  # Need to  extract different information for LSM fits#
  
  param <- lapply(object_list, extract_param, type = type, burnin = burnin)
  param_mcmc_list <- as.mcmc.list(lapply(param, param_to_mcmc))
  
  output <- list(call = match.call())
  if ('psrf' %in% diags) {
    # will omit warning if user omitted the diags argument, and therefore
    # did not explicitly ask for PSRF
    psrf_attrs <- psrf_param(param_mcmc_list, warn_1chain = warn_1chain)
    if (!is.null(psrf_attrs)) {
      output <- c(output, psrf = list(psrf_attrs))
    }
  }
  
  if ('raftery' %in% diags) {
    raft_attrs <- lapply(param_mcmc_list, raftery_param)
    if (!is.null(raft_attrs)) {
      output <- c(output, raftery = list(raft_attrs))
    }
  }
  
  if (('traceplot' %in% diags || 'autocorr' %in% diags) &&
      (nchain(param_mcmc_list) > 1) && missing(col)) {
    chain_dex <- seq_along(param_mcmc_list)
    legend <- paste("Chain", chain_dex, "=", grDevices::palette()[chain_dex],
                    collapse = '\n')
    message("Plot Color Legend:\n", legend)
  }
  
  if ('traceplot' %in% diags) {
    shape_args <- plot_shape(param_mcmc_list, type = type)
    par(mfrow = c(shape_args$nrows, shape_args$ncols), mar = rep(2, 4))
    traceplot(param_mcmc_list[, shape_args$plotdex], 
              col = col, lty = lty)
  }
  
  if ('autocorr' %in% diags) {
    shape_args <- plot_shape(param_mcmc_list, type = type)
    par(mfrow = c(shape_args$nrows, shape_args$ncols), mar = rep(2, 4))
    autocorr_param(param_mcmc_list[, shape_args$plotdex], col = col, lty = lty)
  }
  if (length(output) > 1) {
    class(output) <- "HLSMdiag"
    return(output)
  }
}


# Output Printing ---------------------------------------------------------

call_summary <- function(call) {
  cat("\nCall:\n", paste(deparse(call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
}

print.HLSMdiag <- function(x, ...) {
  # allows this function to be flexible to adding more diagnostic summaries
  summary_funcs <- list(call = call_summary,
                        psrf = psrf_summary,
                        raftery = raftery_summary)
  if (!is(x, 'HLSMdiag')) {
    stop("This function does not work on non-HLSMdiag objects.")
  }
  
  if (!all(names(x) %in% names(summary_funcs))) {
    stop("HLSMdiag function updated without updating print.HLSMdiag.\n",
         "Please contact the maintainer to fix")
  }
  
  for (el in names(x)) {
    func <- summary_funcs[[el]]
    obj <- x[[el]]
    func(obj)
    cat('\n')
  }
  
  cat("Detailed Information:\n")
  cat("To review detailed diagnostic information for each variable,\n")
  cat("access this object as a list with `$`.\n")
}