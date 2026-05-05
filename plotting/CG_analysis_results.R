# CG ERGM Analysis Results
# Written by Carter T. Butts
# Modified by Yuanming Song on September 19, 2025

# Settings - controlled by individual code blocks
maindir <- if (exists("maindir")) maindir else normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
results_dir <- file.path(maindir, "results")
output_dir <- file.path(maindir, "plots")

# Load CG ERGM results
load_cg_results <- function(project_dir) {
  cg_file <- file.path(project_dir, "results", "csh_cssc_1node_fit_CG.Rdata")
  
  if (!file.exists(cg_file)) {
    stop("CG results file not found: ", cg_file)
  }
  
  # Load the data
  env <- new.env()
  load(cg_file, envir = env)
  
  return(list(
    best_fit = env$best.fit,
    best_gof = env$best.gof,
    best_terms = env$best.terms
  ))
}

# Create LaTeX table for CG ERGM results
create_cg_latex_table <- function(best_fit, project_dir, save_file = FALSE) {
  # Extract coefficients and standard errors
  coefs <- best_fit$coef
  ses <- best_fit$se
  
  # Calculate z-values and p-values
  z_values <- coefs / ses
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Create significance codes
  sig_codes <- ifelse(p_values < 0.001, "***",
                     ifelse(p_values < 0.01, "**",
                           ifelse(p_values < 0.05, "*",
                                 ifelse(p_values < 0.1, ".", " "))))
  
  # Get composition distribution from train.sample
  train_sample <- best_fit$train.sample
  comp_dist <- train_sample$type.count
  comp_labels <- sapply(train_sample$types, function(x) {
    if (length(x) == 1) {
      as.character(sum(x))
    } else {
      paste0("c(\"", paste(x, collapse = "\", \""), "\")")
    }
  })
  
  # Start creating LaTeX table
  latex_lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Coarse-Grained ERGM Parameter Estimates}",
    "\\label{tab:cg_ergm}",
    "\\begin{tabular}{lrrrr}",
    "\\hline",
    "Parameter & Estimate & Std.Err & Z value & Pr($>|z|$) \\\\",
    "\\hline"
  )
  
  # Add parameter rows
  param_names <- names(coefs)
  for (i in 1:length(param_names)) {
    param_name <- param_names[i]
    estimate <- sprintf("%.3f", coefs[i])
    std_err <- sprintf("%.3f", ses[i])
    z_val <- sprintf("%.3f", z_values[i])
    
    # Format p-value
    if (p_values[i] < 2.2e-16) {
      p_val <- "$\\textless 2.2 \\times 10^{-16}$"
    } else {
      p_val <- sprintf("%.3e", p_values[i])
    }
    
    latex_lines <- c(latex_lines,
      paste0(param_name, " & ", estimate, " & ", std_err, " & ", z_val, " & ", p_val, " \\\\"))
  }
  
  # Add footer
  latex_lines <- c(latex_lines,
    "\\hline",
    paste0("\\multicolumn{5}{l}{Test deviance ", sprintf("%.1f", best_fit$deviance.test), " (null deviance 0)} \\\\"),
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  
  # Print table to console
  cat("\nCG ERGM Parameter Table:\n")
  cat("========================\n\n")
  
  cat("Composition distribution:\n")
  for (i in 1:length(comp_labels)) {
    cat(sprintf("%-15s %d\n", comp_labels[i], comp_dist[i]))
  }
  
  cat("\nParameter Estimates:\n")
  cat(sprintf("%-25s %8s %8s %8s %10s %s\n", "Parameter", "Estimate", "Std.Err", "Z value", "Pr(>|z|)", ""))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  for (i in 1:length(param_names)) {
    param_name <- param_names[i]
    estimate <- sprintf("%.3f", coefs[i])
    std_err <- sprintf("%.3f", ses[i])
    z_val <- sprintf("%.3f", z_values[i])
    
    if (p_values[i] < 2.2e-16) {
      p_val <- "< 2.2e-16"
    } else {
      p_val <- sprintf("%.3e", p_values[i])
    }
    
    cat(sprintf("%-25s %8s %8s %8s %10s %s\n", param_name, estimate, std_err, z_val, p_val, sig_codes[i]))
  }
  
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat(sprintf("Test deviance %.1f (null deviance 0)\n", best_fit$deviance.test))
  
  # Save to file if requested
  if (save_file) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, "CG_ERGM_table.tex")
    
    # Add LaTeX-ready header comments for coding agent
    final_latex <- c(
      "% CG ERGM Parameter Table",
      "% To include in LaTeX document, use: \\input{CG_ERGM_table.tex}",
      "% This table shows Coarse-Grained ERGM parameter estimates with composition distribution",
      "",
      latex_lines
    )
    
    writeLines(final_latex, output_file)
    cat("\nLaTeX table saved to:", output_file, "\n")
  }
  
  return(latex_lines)
}

# Auto-execute when sourced
project_dir <- maindir

cat("Loading CG ERGM results...\n")
results <- load_cg_results(project_dir)

cat("Creating LaTeX table...\n")
# Use save_plt from Rmd environment (defaults to FALSE if not set)
save_plot <- if(exists("save_plt")) save_plt else FALSE
latex_table <- create_cg_latex_table(results$best_fit, project_dir, save_file = save_plot)

cat("Creating GOF plots...\n")

# Required libraries are loaded in Rmd

# Extract GOF data
degree_obs <- results$best_gof$degree.obs
degree_mean <- results$best_gof$degree.mean
degree_q025 <- results$best_gof$degree.q025
degree_q975 <- results$best_gof$degree.q975

esp_obs <- results$best_gof$esp.obs
esp_mean <- results$best_gof$esp.mean
esp_q025 <- results$best_gof$esp.q025
esp_q975 <- results$best_gof$esp.q975

comp_obs <- results$best_gof$comp.obs
comp_mean <- results$best_gof$comp.mean
comp_q025 <- results$best_gof$comp.q025
comp_q975 <- results$best_gof$comp.q975

# Get system types
type_ind <- results$best_gof$type.ind
types <- results$best_gof$types

# Create system labels mapping
system_labels <- c("0% CSSC", "100% CSSC", "50% CSSC")
type_mapping <- c(1, 2, 3)  # type 1=0%, type 2=100%, type 3=50%

# Function to create GOF subplot for one system
create_gof_subplot <- function(obs_data, mean_data, q025_data, q975_data, 
                              system_type, stat_name, x_labels) {
  
  sys_rows <- which(type_ind == system_type)
  if (length(sys_rows) == 0) return(NULL)
  
  # Get data for this system
  obs_vals <- obs_data[sys_rows, , drop = FALSE]
  mean_vals <- colMeans(mean_data[sys_rows, , drop = FALSE], na.rm = TRUE)
  q025_vals <- colMeans(q025_data[sys_rows, , drop = FALSE], na.rm = TRUE)
  q975_vals <- colMeans(q975_data[sys_rows, , drop = FALSE], na.rm = TRUE)
  obs_means <- colMeans(obs_vals, na.rm = TRUE)
  
  # Create data frame
  plot_data <- data.frame(
    x = 1:length(mean_vals),
    x_labels = x_labels[1:length(mean_vals)],
    sim_mean = mean_vals,
    obs_mean = obs_means,
    q025 = q025_vals,
    q975 = q975_vals
  )
  
  # Add individual observed points
  obs_points <- data.frame()
  for (i in 1:nrow(obs_vals)) {
    obs_row <- data.frame(
      x = 1:ncol(obs_vals),
      y = as.numeric(obs_vals[i, ]),
      sample = i
    )
    obs_points <- rbind(obs_points, obs_row)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x)) +
    # Confidence interval
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "red", alpha = 0.2) +
    # Simulated mean
    geom_line(aes(y = sim_mean), color = "red", linewidth = 1.2) +
    # Observed mean (dashed)
    geom_line(aes(y = obs_mean), color = "black", linetype = "dashed", linewidth = 1.2) +
    # Individual observed points
    geom_point(data = obs_points, aes(x = x, y = y), color = "black", size = 1, alpha = 0.7) +
    scale_x_continuous(breaks = 1:length(mean_vals), labels = x_labels[1:length(mean_vals)]) +
    {if(grepl("Component", stat_name)) scale_y_log10() else NULL} +
    labs(title = system_labels[system_type], 
         x = "", y = if(grepl("Component", stat_name)) "Count (log scale)" else "Vertices") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Function to create legend plot
create_legend_plot <- function() {
  # Create a simple plot for legend
  legend_data <- data.frame(
    x = c(1, 1, 1, 1, 1),
    y = c(5, 4, 3, 2, 1),
    label = c("Observed points", "Observed mean", "Simulated mean", "95% CI", ""),
    color = c("black", "black", "red", "red", "white"),
    linetype = c("solid", "dashed", "solid", "solid", "solid")
  )
  
  p <- ggplot(legend_data, aes(x = x, y = y)) +
    # Observed points
    geom_point(data = legend_data[1,], size = 2, color = "black") +
    geom_text(data = legend_data[1,], aes(label = label), hjust = 0, nudge_x = 0.1, size = 4) +
    # Observed mean (dashed line)
    geom_segment(data = legend_data[2,], aes(x = x-0.05, xend = x+0.05, y = y, yend = y), 
                 color = "black", linetype = "dashed", linewidth = 1.2) +
    geom_text(data = legend_data[2,], aes(label = label), hjust = 0, nudge_x = 0.1, size = 4) +
    # Simulated mean (solid red line)
    geom_segment(data = legend_data[3,], aes(x = x-0.05, xend = x+0.05, y = y, yend = y), 
                 color = "red", linetype = "solid", linewidth = 1.2) +
    geom_text(data = legend_data[3,], aes(label = label), hjust = 0, nudge_x = 0.1, size = 4) +
    # 95% CI (shaded area)
    geom_rect(data = legend_data[4,], aes(xmin = x-0.05, xmax = x+0.05, ymin = y-0.1, ymax = y+0.1), 
              fill = "red", alpha = 0.2, color = NA) +
    geom_text(data = legend_data[4,], aes(label = label), hjust = 0, nudge_x = 0.1, size = 4) +
    xlim(0.8, 2.5) + ylim(0.5, 5.5) +
    theme_void()
  
  return(p)
}

# 1. DEGREE DISTRIBUTION PLOTS (3 systems in one figure)
max_degree <- min(15, ncol(degree_obs) - 1)
degree_labels <- paste0("degree", 0:max_degree)

degree_plots <- list()
for (sys in 1:3) {
  p <- create_gof_subplot(
    degree_obs[, 1:(max_degree+1)], degree_mean[, 1:(max_degree+1)], 
    degree_q025[, 1:(max_degree+1)], degree_q975[, 1:(max_degree+1)],
    sys, "Degree Distribution", degree_labels[1:(max_degree+1)]
  )
  degree_plots[[sys]] <- p
}

# Create legend plot
legend_plot <- create_legend_plot()

degree_combined <- grid.arrange(
  grobs = c(degree_plots, list(legend_plot)), 
  ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  top = "Degree Distribution - Goodness of Fit",
  left = "Vertices",
  bottom = "Degree"
)

# 2. ESP DISTRIBUTION PLOTS (3 systems in one figure)
max_esp <- min(15, ncol(esp_obs) - 1)
esp_labels <- paste0("esp", 0:max_esp)

esp_plots <- list()
for (sys in 1:3) {
  p <- create_gof_subplot(
    esp_obs[, 1:(max_esp+1)], esp_mean[, 1:(max_esp+1)], 
    esp_q025[, 1:(max_esp+1)], esp_q975[, 1:(max_esp+1)],
    sys, "ESP Distribution", esp_labels[1:(max_esp+1)]
  )
  esp_plots[[sys]] <- p
}

esp_combined <- grid.arrange(
  grobs = c(esp_plots, list(legend_plot)), 
  ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  top = "ESP Distribution - Goodness of Fit",
  left = "Vertices", 
  bottom = "Edgewise Shared Partners"
)

# 3. COMPONENT STATISTICS PLOTS (3 systems in one figure)
comp_labels <- colnames(comp_obs)

comp_plots <- list()
for (sys in 1:3) {
  p <- create_gof_subplot(
    comp_obs, comp_mean, comp_q025, comp_q975,
    sys, "Component Statistics", comp_labels
  )
  comp_plots[[sys]] <- p
}

comp_combined <- grid.arrange(
  grobs = c(comp_plots, list(legend_plot)), 
  ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  top = "Component Statistics - Goodness of Fit",
  left = "Count",
  bottom = "Component Type"
)

# Save plots block
save_plots <- if(exists("save_plt")) save_plt else FALSE

if (save_plots) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  png(file.path(output_dir, "CG_degree_GOF.png"), width = 1600, height = 1200, res = 300)
  grid.arrange(
    grobs = c(degree_plots, list(legend_plot)), 
    ncol = 2, nrow = 2,
    layout_matrix = rbind(c(1, 2), c(3, 4)),
    top = "Degree Distribution - Goodness of Fit",
    left = "Vertices",
    bottom = "Degree"
  )
  dev.off()
  
  png(file.path(output_dir, "CG_esp_GOF.png"), width = 1600, height = 1200, res = 300)
  grid.arrange(
    grobs = c(esp_plots, list(legend_plot)), 
    ncol = 2, nrow = 2,
    layout_matrix = rbind(c(1, 2), c(3, 4)),
    top = "ESP Distribution - Goodness of Fit",
    left = "Vertices",
    bottom = "Edgewise Shared Partners"
  )
  dev.off()
  
  png(file.path(output_dir, "CG_comp_GOF.png"), width = 1600, height = 1200, res = 300)
  grid.arrange(
    grobs = c(comp_plots, list(legend_plot)), 
    ncol = 2, nrow = 2,
    layout_matrix = rbind(c(1, 2), c(3, 4)),
    top = "Component Statistics - Goodness of Fit",
    left = "Count",
    bottom = "Component Type"
  )
  dev.off()
  
  cat("CG GOF plots saved to:", output_dir, "\n")
}

cat("CG analysis complete!\n")
