# Simulated CSH-CSSC Network Visualization Script
# Local visualization script for simulated ERGM networks
# SYM, created on Jan 6, 2025

# Required libraries are loaded in Rmd

# Use save_plt from Rmd environment (defaults to FALSE if not set)
save_plot <- if(exists("save_plt")) save_plt else FALSE

# File selection - specify which files to load and plot
# Set to NULL to plot all files, or specify indices like c(1, 3, 5) to plot specific files
files_to_plot <- NULL  # Change this to select specific files, e.g., c(1, 2, 3)

# Simulation type selection
# Set to "both" to plot both AA and CG, "AA" for AA only, "CG" for CG only
simulation_types <- "both"  # Options: "both", "AA", "CG"

# Define text sizes and plot parameters
plot_title_size <- 1
axis_title_size <- 1
axis_text_size <- 1
legend_title_size <- 1
legend_text_size <- 1

# Define local project path
maindir <- if (exists("maindir")) maindir else normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
data_dir <- file.path(maindir, "simulated_networks")
output_dir <- file.path(maindir, "plots")

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to extract CSH and CSSC numbers from filename
extract_node_counts <- function(filename) {
  # Extract pattern like "AA_networks_800CSH_1200CSSC_50sims.Rdata" or "CG_networks_800CSH_1200CSSC_50sims.Rdata"
  aa_pattern <- "AA_networks_(\\d+)CSH_(\\d+)CSSC_\\d+sims\\.Rdata"
  cg_pattern <- "CG_networks_(\\d+)CSH_(\\d+)CSSC_\\d+sims\\.Rdata"
  
  # Try AA pattern first
  aa_matches <- regmatches(filename, regexec(aa_pattern, filename))
  if (length(aa_matches) > 0 && length(aa_matches[[1]]) == 3) {
    n_csh <- as.numeric(aa_matches[[1]][2])
    n_cssc <- as.numeric(aa_matches[[1]][3])
    return(list(n_csh = n_csh, n_cssc = n_cssc, type = "AA"))
  }
  
  # Try CG pattern
  cg_matches <- regmatches(filename, regexec(cg_pattern, filename))
  if (length(cg_matches) > 0 && length(cg_matches[[1]]) == 3) {
    n_csh <- as.numeric(cg_matches[[1]][2])
    n_cssc <- as.numeric(cg_matches[[1]][3])
    return(list(n_csh = n_csh, n_cssc = n_cssc, type = "CG"))
  }
  
  # Fallback: try to extract any numbers from filename
  numbers <- as.numeric(unlist(regmatches(filename, gregexpr("\\d+", filename))))
  if (length(numbers) >= 2) {
    # Determine type from filename
    type <- if (grepl("AA_", filename)) "AA" else if (grepl("CG_", filename)) "CG" else "Unknown"
    return(list(n_csh = numbers[1], n_cssc = numbers[2], type = type))
  } else {
    return(list(n_csh = NA, n_cssc = NA, type = "Unknown"))
  }
}

# Function to plot simulated CSH-CSSC network with colored nodes
plot_simulated_network <- function(net, sim_idx, n_csh, n_cssc, resnames, sim_type = "", filename = "") {
  tryCatch({
    # Get node colors based on Resname
    node_colors <- ifelse(resnames == "CSH", "red", "blue")
  
  # Calculate layout
  # coords <- gplot.layout.fruchtermanreingold(net, NULL)
  coords <- gplot.layout.fruchtermanreingold(net, NULL)  # Force-directed layout
  # coords <- gplot.layout.circle(net, NULL)  # Circle layout for better spacing
  # coords <- gplot.layout.random(net, NULL)
  
  # Create plot with improved spacing
  plot(net, 
       coord = coords,
       vertex.col = node_colors,
       vertex.cex = 0.3,           # Smaller nodes
       vertex.border = "black",
       edge.col = "gray80",         # Lighter edges
       edge.lwd = 0.2,              # Thinner edges
       edge.curve = 0.1,            # Slight curve to reduce overlap
       main = paste0(sim_type, " Network ", sim_idx, " - CSH ", n_csh, ":", n_cssc, " CSSC"),
       cex.main = plot_title_size,
       pad = 0.1)                   # Add padding around plot
  
  # Add legend
  legend("topright", 
         legend = c("CSH", "CSSC"), 
         col = c("red", "blue"), 
         pch = 19, 
         cex = legend_text_size,
         title = "Node Type",
         title.cex = legend_title_size)
  
  # Calculate network statistics
  n_edges <- network.edgecount(net)
  n_nodes <- network.size(net)
  
  # Check if network has edges
  if (n_edges == 0) {
    # Empty network
    largest_comp_pct <- 0
    avg_degree_csh <- 0
    avg_degree_cssc <- 0
    pct_cssc_largest <- 0
  } else {
    # Largest component size as percentage
    components <- component.dist(net)
    largest_comp_size <- max(components$csize)
    largest_comp_pct <- round((largest_comp_size / n_nodes) * 100, 1)
    
    # Average degree for CSH and CSSC/CSX (excluding isolates)
    degrees <- degree(net)
    csh_indices <- which(resnames == "CSH")
    # Handle both CSSC (AA) and CSX (CG) naming
    cssc_indices <- which(resnames == "CSSC" | resnames == "CSX")
    
    if (length(csh_indices) > 0) {
      csh_degrees <- degrees[csh_indices]
      csh_non_isolates <- csh_degrees[csh_degrees > 0]
      if (length(csh_non_isolates) > 0) {
        avg_degree_csh <- round(mean(csh_non_isolates), 2)
      } else {
        avg_degree_csh <- 0
      }
    } else {
      avg_degree_csh <- 0
    }
    
    if (length(cssc_indices) > 0) {
      cssc_degrees <- degrees[cssc_indices]
      cssc_non_isolates <- cssc_degrees[cssc_degrees > 0]
      if (length(cssc_non_isolates) > 0) {
        avg_degree_cssc <- round(mean(cssc_non_isolates), 2)
      } else {
        avg_degree_cssc <- 0
      }
    } else {
      avg_degree_cssc <- 0
    }
    
    # Composition in largest component
    # Find nodes in largest component - use first occurrence if multiple max sizes
    max_comp_size <- max(components$csize)
    largest_comp_id <- which(components$csize == max_comp_size)[1]  # Take first if multiple
    largest_comp_nodes <- which(components$membership == largest_comp_id)
    
    cat("    Largest component ID:", largest_comp_id, "Size:", max_comp_size, "\n")
    cat("    Largest component nodes:", length(largest_comp_nodes), "nodes\n")
    cat("    Node indices:", head(largest_comp_nodes, 5), if(length(largest_comp_nodes) > 5) "..." else "", "\n")
    
    if (length(largest_comp_nodes) > 0) {
      # Get Resname values for largest component nodes
      largest_comp_resnames <- resnames[largest_comp_nodes]
      largest_comp_csh <- sum(largest_comp_resnames == "CSH")
      # Handle both CSSC (AA) and CSX (CG) naming
      largest_comp_cssc <- sum(largest_comp_resnames == "CSSC" | largest_comp_resnames == "CSX")
      largest_comp_total <- largest_comp_csh + largest_comp_cssc
      
      if (largest_comp_total > 0) {
        pct_csh_largest <- round((largest_comp_csh / largest_comp_total) * 100, 1)
        pct_cssc_largest <- round((largest_comp_cssc / largest_comp_total) * 100, 1)
      } else {
        pct_csh_largest <- 0
        pct_cssc_largest <- 0
      }
    } else {
      pct_csh_largest <- 0
      pct_cssc_largest <- 0
    }
  }
  
    # Calculate percent of all CSSC/CSX nodes in largest component
    total_cssc <- sum(resnames == "CSSC" | resnames == "CSX")
    if (total_cssc > 0) {
      pct_cssc_in_largest <- round((largest_comp_cssc / total_cssc) * 100, 1)
    } else {
      pct_cssc_in_largest <- 0
    }
    
    # Create statistics text (each on separate line, positioned at left bottom)
    stats_lines <- c(
      paste0("Largest comp: ", largest_comp_pct, "%"),
      paste0("CSH avg degree: ", avg_degree_csh),
      paste0("CSSC avg degree: ", avg_degree_cssc),
      paste0("Largest comp %CSSC: ", pct_cssc_largest),
      paste0("CSSC in largest: ", pct_cssc_in_largest, "%")
    )
    
    # Add each line at the bottom left
    for (i in 1:length(stats_lines)) {
      mtext(stats_lines[i], side = 1, line = -i, cex = axis_text_size, adj = 0)
    }
    
  }, error = function(e) {
    cat("Error plotting network", sim_idx, ":", e$message, "\n")
    # Create a simple error plot
    plot(1, 1, type = "n", main = paste0("Error in Network ", sim_idx), 
         xlab = "", ylab = "", axes = FALSE)
    text(1, 1, paste("Error:", e$message), cex = 0.8)
  })
}

# Find simulated network files at the beginning
cat("Simulated CSH-CSSC Network Visualization\n")
cat("========================================\n")
cat("Searching for simulated network files...\n")

# Search for both AA and CG files based on simulation_types setting
sim_files <- c()
if (simulation_types == "both" || simulation_types == "AA") {
  aa_files <- list.files(data_dir, pattern = "AA_networks.*\\.Rdata$", full.names = TRUE)
  sim_files <- c(sim_files, aa_files)
}
if (simulation_types == "both" || simulation_types == "CG") {
  cg_files <- list.files(data_dir, pattern = "CG_networks.*\\.Rdata$", full.names = TRUE)
  sim_files <- c(sim_files, cg_files)
}

if (length(sim_files) == 0) {
  cat("No simulated network files found in:", data_dir, "\n")
} else {
  cat("Found", length(sim_files), "simulated network files:\n")
  for (i in 1:length(sim_files)) {
    cat("  ", i, ":", basename(sim_files[i]), "\n")
  }
  
  # Select files to plot
  if (is.null(files_to_plot)) {
    files_to_plot <- 1:length(sim_files)
    cat("Plotting all files\n")
  } else {
    cat("Plotting files:", paste(files_to_plot, collapse = ", "), "\n")
  }
  
  # Process each selected file
  for (file_idx in files_to_plot) {
    if (file_idx > length(sim_files)) {
      cat("Warning: File index", file_idx, "out of range\n")
      next
    }
    
    sim_file <- sim_files[file_idx]
    cat("\nLoading:", basename(sim_file), "\n")
    
    # Extract node counts and type from filename
    node_counts <- extract_node_counts(basename(sim_file))
    filename_csh <- node_counts$n_csh
    filename_cssc <- node_counts$n_cssc
    sim_type <- node_counts$type
    
    # Load the simulated networks
    load(sim_file)
    
    # Check if simulated_networks exists
    if (!exists("simulated_networks")) {
      cat("Error: simulated_networks object not found in file\n")
      next
    }
    
    n_networks <- length(simulated_networks)
    cat("Found", n_networks, "simulated networks\n")
    
    # Visualize only the first network from this file
    n_to_plot <- min(1, n_networks)
    cat("Visualizing first", n_to_plot, "network from this file...\n")
    
    for (i in 1:n_to_plot) {
      cat("  Plotting network", i, "of", n_to_plot, "...\n")
      
      tryCatch({
        net <- simulated_networks[[i]]
        
        # Extract Resname information once per network
        resnames <- sapply(net$val, function(x) x$Resname)
        
        # Use filename counts for title, actual counts for analysis
        actual_csh <- sum(resnames == "CSH")
        actual_cssc <- sum(resnames == "CSSC")
        
        # Use filename counts for title if available, otherwise use actual counts
        title_csh <- if (!is.na(filename_csh)) filename_csh else actual_csh
        title_cssc <- if (!is.na(filename_cssc)) filename_cssc else actual_cssc
      
      if (save_plot) {
        # Create filename with simulation type (single network per type)
        filename <- paste0("simulated_network_", sim_type, "_fitted.png")
        filepath <- file.path(output_dir, filename)
        
        # Save plot - twice as wide
        png(filepath, width = 3200, height = 1600, res = 400)
        par(mar = c(1, 1, 2, 1))
        
        plot_simulated_network(net, i, title_csh, title_cssc, resnames, sim_type, basename(sim_file))
        
        dev.off()
        
        cat("    Plot saved:", filepath, "\n")
      } else {
        # Just display the plot
        par(mar = c(1, 1, 2, 1))
        plot_simulated_network(net, i, title_csh, title_cssc, resnames, sim_type, basename(sim_file))
      }
      
      # Calculate and display the same statistics with error handling
      n_edges <- network.edgecount(net)
      n_nodes <- network.size(net)
      
      if (n_edges == 0) {
        # Empty network
        largest_comp_pct <- 0
        avg_degree_csh <- 0
        avg_degree_cssc <- 0
        pct_cssc_largest <- 0
      } else {
        components <- component.dist(net)
        largest_comp_size <- max(components$csize)
        largest_comp_pct <- round((largest_comp_size / n_nodes) * 100, 1)
        
        degrees <- degree(net)
        csh_indices <- which(resnames == "CSH")
        # Handle both CSSC (AA) and CSX (CG) naming
        cssc_indices <- which(resnames == "CSSC" | resnames == "CSX")
        
        if (length(csh_indices) > 0) {
          csh_degrees <- degrees[csh_indices]
          csh_non_isolates <- csh_degrees[csh_degrees > 0]
          if (length(csh_non_isolates) > 0) {
            avg_degree_csh <- round(mean(csh_non_isolates), 2)
          } else {
            avg_degree_csh <- 0
          }
        } else {
          avg_degree_csh <- 0
        }
        
        if (length(cssc_indices) > 0) {
          cssc_degrees <- degrees[cssc_indices]
          cssc_non_isolates <- cssc_degrees[cssc_degrees > 0]
          if (length(cssc_non_isolates) > 0) {
            avg_degree_cssc <- round(mean(cssc_non_isolates), 2)
          } else {
            avg_degree_cssc <- 0
          }
        } else {
          avg_degree_cssc <- 0
        }
        
        # Find largest component - use first occurrence if multiple max sizes
        max_comp_size <- max(components$csize)
        largest_comp_id <- which(components$csize == max_comp_size)[1]  # Take first if multiple
        largest_comp_nodes <- which(components$membership == largest_comp_id)
        
        if (length(largest_comp_nodes) > 0) {
          # Get Resname values for largest component nodes
          largest_comp_resnames <- resnames[largest_comp_nodes]
          largest_comp_csh <- sum(largest_comp_resnames == "CSH")
          # Handle both CSSC (AA) and CSX (CG) naming
          largest_comp_cssc <- sum(largest_comp_resnames == "CSSC" | largest_comp_resnames == "CSX")
          largest_comp_total <- largest_comp_csh + largest_comp_cssc
          
          if (largest_comp_total > 0) {
            pct_cssc_largest <- round((largest_comp_cssc / largest_comp_total) * 100, 1)
          } else {
            pct_cssc_largest <- 0
          }
        } else {
          pct_cssc_largest <- 0
        }
      }
      
        # Calculate percent of all CSSC/CSX nodes in largest component
        total_cssc <- sum(resnames == "CSSC" | resnames == "CSX")
        if (total_cssc > 0) {
          pct_cssc_in_largest <- round((largest_comp_cssc / total_cssc) * 100, 1)
        } else {
          pct_cssc_in_largest <- 0
        }
        
        cat("    Largest comp:", largest_comp_pct, "% | CSH avg degree:", avg_degree_csh, 
            "| CSSC avg degree:", avg_degree_cssc, "| Largest comp %CSSC:", pct_cssc_largest, 
            "| CSSC in largest:", pct_cssc_in_largest, "%\n")
        
      }, error = function(e) {
        cat("    Error processing network", i, ":", e$message, "\n")
      })
    }
    
    cat("Completed file:", basename(sim_file), "\n")
  }
  
  cat("\nSimulated network visualization complete!\n")
  if (save_plot) {
    cat("Plots saved in:", output_dir, "\n")
  } else {
    cat("Plots displayed in R (not saved)\n")
  }
}

cat("\nVisualization complete!\n")
