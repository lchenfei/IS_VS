#' =============================================================================
#' Main Results Runner - IJDS Paper
#' =============================================================================
#' 
#' Single script to generate the main result table across all 4 examples.
#' Runs the 7 core methods: IS-VS, IS-CE, IS-Pareto, WAMK-SIS, Lasso, SpAM, RF-RFE
#' 
#' Note: Lasso is not available for Example 3 (D=10), so 6 methods run for Ex3.
#' 
#' Usage:
#'   source("run_main_results.R")
#'   
#'   # Run all examples with all methods
#'   results <- run_all_experiments()
#'   
#'   # Run specific example
#'   results <- run_example(3, n_reps = 25, ncores = 25)
#'   
#'   # Run specific method on specific example
#'   results <- run_single(example = 1, method = "IS-VS", n_reps = 25, ncores = 25)
#'   
#'   # Generate summary table
#'   summary_table <- generate_summary_table(results)
#' 
#' =============================================================================

library(doParallel)

# =============================================================================
# Configuration
# =============================================================================

ROOT <- getwd()
RESULTS_DIR <- file.path(ROOT, "results")

if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Core methods for main results table (7 methods, but Lasso excluded for Ex3)
CORE_METHODS <- c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")

# Example configurations
EXAMPLE_CONFIG <- list(
  Ex1 = list(
    D = 4,
    l_target = 18.99,
    n_pairs = 6,
    seed_base = 123,
    sampling = "independent",
    functions_file = "functions_ex1.R",
    methods = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")
  ),
  Ex2 = list(
    D = 5,
    l_target = 24.2,
    n_pairs = 10,
    seed_base = 1235,
    sampling = "correlated",
    functions_file = "functions_ex2.R",
    methods = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")
  ),
  Ex3 = list(
    D = 10,
    l_target = 18.74,
    n_pairs = 45,
    seed_base = 123,
    sampling = "independent",
    functions_file = "functions_ex3.R",
    methods = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")  
  ),
  Ex4 = list(
    D = 50,
    l_target = 18.42,
    n_pairs = 1225,
    seed_base = 127,
    sampling = "independent",
    functions_file = "functions_ex4.R",
    methods = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")
  )
)

# Method file mapping
METHOD_FILES <- list(
  "IS-VS"     = "IS_VS.R",
  "IS-CE"     = "IS_CE.R",
  "IS-Pareto" = "IS_Pareto.R",
  "WAMK-SIS"  = "wamk.R",
  "Lasso"     = "Lasso.R",
  "SpAM"      = "SpAM.R",
  "RF-RFE"    = "randomForest.R",
  "OptiTree" = "OptiTreeStrata.R" 
)

METHOD_FUNCTIONS <- list(
  "IS-VS"     = "run_is_vs",
  "IS-CE"     = "run_is_ce",
  "IS-Pareto" = "run_is_pareto",
  "WAMK-SIS"  = "run_wamk_sis",
  "Lasso"     = "run_lasso",
  "SpAM"      = "run_spam",
  "RF-RFE"    = "run_rf_rfe", 
  "OptiTree" = "run_optitree"
)

# =============================================================================
# Core Functions
# =============================================================================

#' Run a single method on a single example
#' 
#' @param example Integer 1-4 specifying which example
#' @param method Character string specifying method name
#' @param n_reps Number of replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param save_results Whether to save results to file (default: TRUE)
#' @return List with results
run_single <- function(example, 
                       method, 
                       n_reps = 25, 
                       ncores = 25,
                       save_results = TRUE) {
  
  ex_name <- paste0("Ex", example)
  config <- EXAMPLE_CONFIG[[ex_name]]
  
  # Check if method is available for this example
 if (!method %in% config$methods) {
    warning(sprintf("Method '%s' not available for %s. Skipping.", method, ex_name))
    return(NULL)
  }
  
  message(sprintf("\n=== Running %s on %s ===", method, ex_name))
  message(sprintf("D=%d, l_target=%.2f, n_reps=%d, ncores=%d", 
                  config$D, config$l_target, n_reps, ncores))
  
  # Source the method file
  method_file <- file.path(ROOT, ex_name, "methods", 
                           paste0("ex", example, "_", METHOD_FILES[[method]]))
  
  if (!file.exists(method_file)) {
    stop("Method file not found: ", method_file)
  }
  
  source(method_file)
  
  # Get and run the function
  run_fn <- get(METHOD_FUNCTIONS[[method]])
  
  start_time <- Sys.time()
  results <- run_fn(n_reps = n_reps, ncores = ncores, ROOT = ROOT)
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  
  message(sprintf("Completed in %s", format(elapsed)))
  
  # Extract estimates
  if (is.list(results) && !is.null(results[[1]]$est_vec)) {
    # Results are in list format with est_vec
    estimates <- sapply(results, function(x) mean(x$est_vec))
  } else if (is.list(results)) {
    # Results might be raw POE_est format
    estimates <- unlist(results)
  } else {
    estimates <- results
  }
  
  # Package output
  output <- list(
    example = ex_name,
    method = method,
    estimates = estimates,
    mean_est = mean(estimates, na.rm = TRUE),
    std_est = sd(estimates, na.rm = TRUE),
    n_reps = length(estimates),
    elapsed_time = elapsed,
    config = config,
    timestamp = Sys.time()
  )
  
  # Save results
  if (save_results) {
    output_file <- file.path(RESULTS_DIR, 
                             sprintf("%s_%s.Rdata", ex_name, gsub("-", "_", method)))
    save(output, file = output_file)
    message("Results saved to: ", output_file)
  }
  
  return(output)
}

#' Run all methods for a single example
#' 
#' @param example Integer 1-4 specifying which example
#' @param methods Vector of method names (default: all available for example)
#' @param n_reps Number of replications
#' @param ncores Number of parallel cores
#' @param save_results Whether to save individual results
#' @return List of results for each method
run_example <- function(example, 
                        methods = NULL, 
                        n_reps = 25, 
                        ncores = 25,
                        save_results = TRUE) {
  
  ex_name <- paste0("Ex", example)
  config <- EXAMPLE_CONFIG[[ex_name]]
  
  if (is.null(methods)) {
    methods <- config$methods
  }
  
  message(sprintf("\n########## Running %s: %d methods ##########", ex_name, length(methods)))
  
  results <- list()
  
  for (method in methods) {
    tryCatch({
      results[[method]] <- run_single(example, method, n_reps, ncores, save_results)
    }, error = function(e) {
      warning(sprintf("Method %s failed on %s: %s", method, ex_name, e$message))
      results[[method]] <<- list(
        example = ex_name,
        method = method,
        error = e$message
      )
    })
  }
  
  return(results)
}

#' Run all experiments (all methods on all examples)
#' 
#' @param examples Vector of example numbers (default: 1:4)
#' @param methods Vector of method names (default: all core methods)
#' @param n_reps Number of replications
#' @param ncores Number of parallel cores
#' @param save_results Whether to save results
#' @return Nested list of all results
run_all_experiments <- function(examples = 1:4,
                                methods = NULL,
                                n_reps = 25,
                                ncores = 25,
                                save_results = TRUE) {
  
  message("=============================================================")
  message("         IJDS Paper - Main Results Generation")
  message("=============================================================")
  message(sprintf("Examples: %s", paste(examples, collapse = ", ")))
  message(sprintf("Replications: %d, Cores: %d", n_reps, ncores))
  message("=============================================================\n")
  
  all_start <- Sys.time()
  all_results <- list()
  
  for (ex in examples) {
    ex_name <- paste0("Ex", ex)
    ex_methods <- methods
    if (is.null(ex_methods)) {
      ex_methods <- EXAMPLE_CONFIG[[ex_name]]$methods
    } else {
      # Filter to only available methods for this example
      ex_methods <- intersect(ex_methods, EXAMPLE_CONFIG[[ex_name]]$methods)
    }
    
    all_results[[ex_name]] <- run_example(ex, ex_methods, n_reps, ncores, save_results)
  }
  
  all_end <- Sys.time()
  total_time <- all_end - all_start
  
  message("\n=============================================================")
  message(sprintf("All experiments completed in %s", format(total_time)))
  message("=============================================================")
  
  # Save combined results
  if (save_results) {
    combined_file <- file.path(RESULTS_DIR, "all_results.Rdata")
    save(all_results, file = combined_file)
    message("Combined results saved to: ", combined_file)
  }
  
  return(all_results)
}

#' Generate summary table from results
#' 
#' @param results Results from run_all_experiments or run_example
#' @param digits Number of decimal places for formatting
#' @return Data frame with summary statistics
generate_summary_table <- function(results, digits = 6) {
  
  # Determine structure of results
  if (!is.null(results$Ex1) || !is.null(results$Ex2)) {
    # Results from run_all_experiments (nested by example)
    examples <- names(results)
  } else {
    # Results from single example
    examples <- NULL
    results <- list(Single = results)
  }
  
  rows <- list()
  
  for (ex_name in names(results)) {
    ex_results <- results[[ex_name]]
    
    for (method_name in names(ex_results)) {
      res <- ex_results[[method_name]]
      
      if (is.null(res) || !is.null(res$error)) {
        next
      }
      
      rows[[length(rows) + 1]] <- data.frame(
        Example = res$example,
        Method = res$method,
        Mean = round(res$mean_est, digits),
        SD = round(res$std_est, digits),
        N = res$n_reps,
        Time = format(res$elapsed_time),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0) {
    return(data.frame())
  }
  
  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL
  
  return(summary_df)
}

#' Print formatted results table (LaTeX-ready)
#' 
#' @param results Results from run_all_experiments
#' @param format Output format: "console", "latex", or "markdown"
print_results_table <- function(results, format = "console") {
  
  summary_df <- generate_summary_table(results)
  
  if (nrow(summary_df) == 0) {
    message("No results to display.")
    return(invisible(NULL))
  }
  
  if (format == "latex") {
    # Create LaTeX table
    cat("\\begin{table}[htbp]\n")
    cat("\\centering\n")
    cat("\\caption{Probability of Exceedance Estimates}\n")
    cat("\\begin{tabular}{llcccc}\n")
    cat("\\hline\n")
    cat("Example & Method & Mean & SD & N & Time \\\\\n")
    cat("\\hline\n")
    
    for (i in 1:nrow(summary_df)) {
      row <- summary_df[i, ]
      cat(sprintf("%s & %s & %.6f & %.6f & %d & %s \\\\\n",
                  row$Example, row$Method, row$Mean, row$SD, row$N, row$Time))
    }
    
    cat("\\hline\n")
    cat("\\end{tabular}\n")
    cat("\\end{table}\n")
    
  } else if (format == "markdown") {
    # Create Markdown table
    cat("| Example | Method | Mean | SD | N | Time |\n")
    cat("|---------|--------|------|-----|---|------|\n")
    
    for (i in 1:nrow(summary_df)) {
      row <- summary_df[i, ]
      cat(sprintf("| %s | %s | %.6f | %.6f | %d | %s |\n",
                  row$Example, row$Method, row$Mean, row$SD, row$N, row$Time))
    }
    
  } else {
    # Console format
    print(summary_df)
  }
  
  return(invisible(summary_df))
}

#' Create wide-format table (methods as columns, examples as rows)
#' 
#' @param results Results from run_all_experiments
#' @param metric Which metric to display: "mean", "sd", or "both"
#' @return Wide-format data frame
create_wide_table <- function(results, metric = "both") {
  
  summary_df <- generate_summary_table(results)
  
  if (nrow(summary_df) == 0) {
    return(data.frame())
  }
  
  # Get unique examples and methods in order
  examples <- c("Ex1", "Ex2", "Ex3", "Ex4")
  examples <- examples[examples %in% unique(summary_df$Example)]
  
  methods <- c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")
  methods <- methods[methods %in% unique(summary_df$Method)]
  
  # Create wide table
  wide_df <- data.frame(Example = examples, stringsAsFactors = FALSE)
  
  for (method in methods) {
    col_values <- c()
    for (ex in examples) {
      row <- summary_df[summary_df$Example == ex & summary_df$Method == method, ]
      if (nrow(row) == 0) {
        col_values <- c(col_values, NA)
      } else if (metric == "mean") {
        col_values <- c(col_values, row$Mean)
      } else if (metric == "sd") {
        col_values <- c(col_values, row$SD)
      } else {
        col_values <- c(col_values, sprintf("%.4f (%.4f)", row$Mean, row$SD))
      }
    }
    wide_df[[method]] <- col_values
  }
  
  return(wide_df)
}

# =============================================================================
# Quick Start Examples
# =============================================================================

message("\n=== Main Results Runner Loaded ===")
message("Quick start examples:")
message("  run_single(1, 'IS-VS')           # Run IS-VS on Example 1")
message("  run_example(2)                   # Run all methods on Example 2")
message("  run_all_experiments()            # Run everything")
message("  generate_summary_table(results)  # Create summary table")
message("  print_results_table(results, 'latex')  # LaTeX output")
message("  create_wide_table(results)       # Wide format for paper")
