# CSH-CSSC Network Visualization Script
# Local visualization script for CSH-CSSC molecular dynamics network data
# Modified from peptide visualization reference
# SYM, created on Jan 6, 2025

# Load required libraries
library(sna)
library(network)

# Settings
save_plot <- FALSE  # Set to TRUE to save plots, FALSE to just display

# Define text sizes and plot parameters
plot_title_size <- 12
axis_title_size <- 10
axis_text_size <- 8
legend_title_size <- 10
legend_text_size <- 8

# Define local project path
maindir <- if (exists("maindir")) maindir else normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
data_dir <- file.path(maindir, "data", "CSH-CSSC_simulation")

# Function to create network from edgelist with node attributes
create_csh_cssc_network <- function(node_data, edgelist, frame_idx) {
  # Get the specific frame's edgelist
  el <- edgelist[[frame_idx]]
  
  # Get number of nodes from node data
  n_nodes <- nrow(node_data)
  
  # Create network object
  net <- network.initialize(n_nodes, directed = FALSE)
  
  # Add edges for this frame
  if (nrow(el) > 0) {
    for (i in 1:nrow(el)) {
      net[el[i,1], el[i,2]] <- 1
    }
  }
  
  # Add vertex attributes
  net %v% "Index" <- node_data$Index
  net %v% "Resname" <- node_data$Resname
  net %v% "Resid" <- node_data$Resid
  
  return(net)
}

# Function to plot CSH-CSSC network with colored nodes
plot_csh_cssc_network <- function(net, composition, system_type, frame_idx, time_ns) {
  # Get node colors based on Resname
  # Note: CG system uses "CSX" but should be treated as "CSSC" for visualization
  resnames <- net %v% "Resname"
  node_colors <- ifelse(resnames == "CSH", "red", 
                       ifelse(resnames == "CSSC" | resnames == "CSX", "blue", "gray"))
  
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
    
    # Average degree for CSH and CSSC (excluding isolates)
    degrees <- degree(net)
    csh_indices <- which(resnames == "CSH")
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
    
    if (length(largest_comp_nodes) > 0) {
      # Get Resname values for largest component nodes
      largest_comp_resnames <- resnames[largest_comp_nodes]
      largest_comp_csh <- sum(largest_comp_resnames == "CSH")
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
  
  # Create layout using Fruchterman-Reingold
  layout_coords <- gplot.layout.fruchtermanreingold(
    net,
    layout.par = list(
      area = 10,
      repulse = 8000,
      niter = 5000
    )
  )
  
  # Calculate actual node counts for title
  actual_csh <- sum(resnames == "CSH")
  actual_cssc <- sum(resnames == "CSSC" | resnames == "CSX")
  
  # Create title with node counts
  title_text <- paste0(system_type, " - ", composition, "% CSSC - CSH ", actual_csh, ":", actual_cssc, " CSSC")
  
  # Plot network
  gplot(net,
        main = title_text,
        cex.main = plot_title_size/10,
        vertex.cex = 1.2,
        vertex.col = node_colors,
        edge.col = "darkgray",
        edge.lwd = 1.5,
        mode = "fruchtermanreingold",
        coord = layout_coords,
        displaylabels = FALSE,
        displayisolates = TRUE,
        usearrows = FALSE)
  
  # Add legend
  legend("topright", 
         legend = c("CSH", "CSSC"),
         col = c("red", "blue"),
         pch = 19,
         cex = legend_text_size/10,
         title = "Node Type",
         title.cex = legend_title_size/10)
  
  # Calculate percent of all CSSC nodes in largest component
  total_cssc <- sum(resnames == "CSSC" | resnames == "CSX")
  if (total_cssc > 0) {
    pct_cssc_in_largest <- round((largest_comp_cssc / total_cssc) * 100, 1)
  } else {
    pct_cssc_in_largest <- 0
  }
  
  # Add statistics text (each on separate line, positioned at left bottom)
  stats_lines <- c(
    paste0("Largest comp: ", largest_comp_pct, "%"),
    paste0("CSH avg degree: ", avg_degree_csh),
    paste0("CSSC avg degree: ", avg_degree_cssc),
    paste0("Largest comp %CSSC: ", pct_cssc_largest),
    paste0("CSSC in largest: ", pct_cssc_in_largest, "%")
  )
  
  # Add each line at the bottom left
  for (i in 1:length(stats_lines)) {
    mtext(stats_lines[i], side = 1, line = -i, cex = axis_text_size/10, adj = 0)
  }
}

# Function to get data paths for each system and composition
get_csh_cssc_paths <- function(system_type, composition) {
  if (system_type == "Atomistic") {
    base_path <- file.path(data_dir, "Fully_Atomistic_32_CSH_equivalents")
  } else if (system_type == "CG") {
    base_path <- file.path(data_dir, "CG_321_CSH_equivalents")
  }
  
  comp_path <- file.path(base_path, paste0(composition, "%"), "1node")
  
  # Get file names based on composition
  if (composition == 0) {
    if (system_type == "Atomistic") {
      nodes_file <- "nodes_CSH32_50mM_1node.rda"
      edges_file <- "equali_csh32_50mM_every10_1node.edgel.stack.rda"
    } else {
      nodes_file <- "nodes_CSH_50mM_cg_1node.rda"
      edges_file <- "equali_csh_50mM_cg_1node.edgel.stack.rda"
    }
  } else if (composition == 25) {
    nodes_file <- "nodes_CSH_CSSC_3to1_50mM_1node.rda"
    edges_file <- "equali_csh_cssc_3to1_50mM_every10_1node.edgel.stack.rda"
  } else if (composition == 50) {
    if (system_type == "Atomistic") {
      nodes_file <- "nodes_CSH_CSSC_1to1_50mM_1node.rda"
      edges_file <- "equali_csh_cssc_1to1_50mM_every10_1node.edgel.stack.rda"
    } else {
      nodes_file <- "nodes_CSH_CSSC_1to1_50mM_cg_1node.rda"
      edges_file <- "equali_csh_cssc_1to1_50mM_cg_1node.edgel.stack.rda"
    }
  } else if (composition == 75) {
    if (system_type == "Atomistic") {
      nodes_file <- "nodes_CSH_CSSC_1to3_50mM_1node_1node.rda"
    } else {
      nodes_file <- "nodes_CSH_CSSC_1to3_50mM_1node.rda"  # CG doesn't have the duplicate
    }
    edges_file <- "equali_csh_cssc_1to3_50mM_every10_1node.edgel.stack.rda"
  } else if (composition == 100) {
    if (system_type == "Atomistic") {
      nodes_file <- "nodes_CSSC_50mM_1node.rda"
      edges_file <- "equali_cssc16_50mM_every10_1node.edgel.stack.rda"
    } else {
      nodes_file <- "nodes_CSSC_50mM_cg_1node.rda"
      edges_file <- "equali_cssc_50mM_cg_1node.edgel.stack.rda"
    }
  }
  
  return(list(
    nodes = file.path(comp_path, nodes_file),
    edges = file.path(comp_path, edges_file)
  ))
}

# Function to load and plot one random network for each system/composition
plot_random_networks <- function() {
  # Define systems and compositions
  systems <- c("Atomistic", "CG")
  atomistic_comps <- c(0, 25, 50, 75, 100)
  cg_comps <- c(0, 50, 100)  # CG only has these compositions
  
  # Create output directory for manuscript figures
  output_dir <- file.path(maindir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set random seed for reproducibility
  set.seed(42)
  
  for (system in systems) {
    compositions <- if (system == "Atomistic") atomistic_comps else cg_comps
    
    cat("\nProcessing", system, "system...\n")
    
    for (comp in compositions) {
      cat("  Loading", comp, "% CSSC composition...\n")
      
      # Get file paths
      paths <- get_csh_cssc_paths(system, comp)
      
      # Check if files exist
      if (!file.exists(paths$nodes) || !file.exists(paths$edges)) {
        cat("    Warning: Files not found for", system, comp, "%\n")
        cat("    Nodes:", paths$nodes, "\n")
        cat("    Edges:", paths$edges, "\n")
        next
      }
      
      # Load data
      load(paths$nodes)  # This loads the node data frame
      load(paths$edges)  # This loads the edgelist
      
      # Get the actual variable names (they vary by file)
      node_vars <- ls(pattern = "nodes_.*")
      edge_vars <- ls(pattern = ".*_.*_.*")
      edge_vars <- edge_vars[!grepl("nodes_", edge_vars)]
      
      if (length(node_vars) == 0 || length(edge_vars) == 0) {
        cat("    Warning: Could not find data variables\n")
        next
      }
      
      # Get the data
      node_data <- get(node_vars[1])
      edgelist <- get(edge_vars[1])
      
      # Randomly select a frame (excluding first and last few frames)
      n_frames <- length(edgelist)
      if (n_frames < 10) {
        frame_idx <- sample(1:n_frames, 1)
      } else {
        frame_idx <- sample(5:(n_frames-5), 1)  # Avoid early/late frames
      }
      
      # Calculate time in nanoseconds (assuming 0.25 ns per frame as in peptide script)
      time_ns <- round(0.25 * (frame_idx - 1), 2)
      
      # Create network
      net <- create_csh_cssc_network(node_data, edgelist, frame_idx)
      
      if (save_plot) {
        # Create filename
        filename <- paste0("network_", system, "_", comp, "pct_CSSC.png")
        filepath <- file.path(output_dir, filename)
        
        # Save plot
        png(filepath, width = 1600, height = 1600, res = 400)
        par(mar = c(1, 1, 2, 1))
        
        plot_csh_cssc_network(net, comp, system, frame_idx, time_ns)
        
        dev.off()
        
        cat("    Plot saved:", filepath, "\n")
      } else {
        # Just display the plot
        par(mar = c(1, 1, 2, 1))
        plot_csh_cssc_network(net, comp, system, frame_idx, time_ns)
      }
      
      cat("    Frame", frame_idx, "of", n_frames, "- Edges:", network.edgecount(net), "\n")
      
      # Clean up loaded variables
      rm(list = c(node_vars, edge_vars))
    }
  }
  
  cat("\nAll network plots completed!\n")
  if (save_plot) {
    cat("Plots saved in:", output_dir, "\n")
  } else {
    cat("Plots displayed in R (not saved)\n")
  }
}

# Main execution
cat("CSH-CSSC Network Visualization\n")
cat("==============================\n")
cat("Loading and plotting random networks from each system and composition...\n")

plot_random_networks()

cat("\nVisualization complete!\n")
