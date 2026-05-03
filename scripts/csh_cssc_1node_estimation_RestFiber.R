#Script to fit network Hamiltonian models to the CSH/CSSC restrained fiber system.
#This version processes restrained fiber simulation data with different structure:
#- Mix system: first 411 nodes are CSSC, rest are CSH
#- Pure system: all CSH nodes
#- Data comes as edgelists in .rda files rather than processed network objects
#Forward selection by held-out deviance is used to choose the model, with parameter
#estimation using a multi-compositional version of Yin-Butts pooling.
#
#Based on CG script, modified for restrained fiber data by SYM on 2025-09-25
#
#To perform estimation, it should suffice to adjust the settings in the modifiable
#arguments section (if desired), and run the code below it.  Some R packages are
#requred; see the libraries and imports section for what is loaded.  See the
#process_MD_data.R script for how to




# get your data into the requisite format for
#this script.


#==============================================================================
#CONFIGURATION PARAMETERS - MODIFY THESE AS NEEDED
#==============================================================================

#Random number seed
base.seed<-456885505

#Base directory and file paths
maindir <- Sys.getenv("NHM_MAIN_DIR", unset = normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE))
base_dir <- maindir

#Directories in which stuff should be put
datdir<-Sys.getenv("NHM_DATA_DIR", unset = file.path(base_dir, "data", "Restrained_Fiber_simulation"))    #Location of restrained fiber data files
outdir<-file.path(base_dir, "results")      #Place to put results
logdir<-file.path(base_dir, "results")      #Place to put log files

#Input data files for restrained fiber simulation (using long trajectories only)
#Mix system file (CSSC + CSH) - longer trajectory with more frames
mix_file<-file.path(datdir, "mix_27A_fibre_edgelist_1to43.rda")
#Pure system file (CSSC only) - longer trajectory with more frames
pure_file<-file.path(datdir, "pure_27A_fibre_edgelist_1to33.rda")

#Output model file, to which fitted model information should be written.
modfile<-file.path(outdir, "csh_cssc_1node_fit_RestFiber.Rdata")

#Temporary checkpoint file for failsafe recovery
checkpoint_file<-file.path(outdir, "csh_cssc_1node_fit_RestFiber_checkpoint.Rdata")

#Intermediate results file for periodic saves (during model selection)
intermediate_file <- file.path(outdir, "csh_cssc_1node_fit_RestFiber_intermediate.Rdata")

#Batch checkpoint file for parallel processing recovery
batch_checkpoint_file <- file.path(outdir, "csh_cssc_1node_fit_RestFiber_batch_checkpoint.Rdata")

#Log file for detailed output
logfile<-file.path(logdir, "csh_cssc_1node_fit_RestFiber.log")

#Log directory for batch logs
batch_log_dir <- file.path(logdir, "RestFiberLog")

#Sample size ratios (calculated from original AA script: 500 train, 250 test, 100 gof = 850 total)
train_ratio <- 500/850  #Training set ratio (58.8% of available data)
test_ratio <- 250/850   #Test set ratio (29.4% of available data) 
gof_ratio <- 100/850    #GOF set ratio (11.8% of available data)
#Note: Ratios sum to 100% and match original AA proportions

#Sample sizes will be calculated dynamically after data loading based on:
#- Minimum available networks across all conditions
#- Above ratios applied to ensure balanced design
#- Conservative buffer to prevent sampling errors

#Estimation parameters for the pooled stochastic approximation algorithm
thin<-1e4               #Thinning interval for MCMC draws; want high enough for decent mixing
subphases<-6            #Number of subphases to use during SA; more=more refined, but much slower
a1<-0.1                 #SA learning rate; too low=too slow, too high=unstable
se.accept<-2            #Number of SEs improvement needed to accept a new model
se.safety<-20           #Safety margin for deviance SE; if the estimated SE of the estimated
                        #deviance is larger than this many multiples of the baseline model,
                        #the candidate model is rejected out of hand (because the MCMC is
                        #not mixing, the model is probably degenerate, and the deviance 
                        #estimate is probably biased)
cores<-10                  #Number of cores to use during estimation; process candidates in batches of 6

#Periodic saving parameters
save_interval <- 2      #Save intermediate results every N terms tested (for monitoring progress)

#Candidate NHM terms.  These may use the "Resname" attribute (and the resnames vector); each must
#be a single term known to ergm or one of the other loaded libraries (see the LIBRARIES) section
#below.  Note that a baseline edge term (and a Krivitsky offset) will be automatically included,
#and should not be listed here.  Be aware that any terms should be well-defined across all 
#compositions (so e.g. you don't want something that can't be defined if all of the nodes are of
#the same type).
resnames<-c("CSH","CSSC")  # Residue types (CSSC nodes in CG systems)
cand<-c(
  'kstar(2)',
  'nodematch("Resname")',
  'nodematch("Resname",levels=I(resnames),diff=TRUE)',
  'nodefactor("Resname",levels=I(resnames)[2])',
  'isolates',
  'dimers',
  'degree(0,by="Resname",homophily=FALSE,levels=I(resnames)[1])',
  'degree(0,by="Resname",homophily=FALSE,levels=I(resnames)[2])',
  'degree(1,by="Resname",homophily=FALSE,levels=I(resnames)[1])',
  'degree(1,by="Resname",homophily=FALSE,levels=I(resnames)[2])',
  'degree(2,by="Resname",homophily=FALSE,levels=I(resnames)[1])',  # Higher degree terms for CSH to broaden distribution
  'degree(2,by="Resname",homophily=FALSE,levels=I(resnames)[2])',  # Higher degree terms for CSSC to broaden distribution
  'degree(3,by="Resname",homophily=FALSE,levels=I(resnames)[2])',  # Even higher degree for CSSC to counter narrow distribution
  'gwdegree(0.25,fixed=TRUE)',
  'gwdegree(0.5,fixed=TRUE)',
  'gwdegree(1,fixed=TRUE)',
  'gwdegree(2,fixed=TRUE)',
  'gwdegree(3,fixed=TRUE)',  # Higher gwdegree parameter to allow more degree heterogeneity
  'components',
  'compsizesum(pow=2)',
  'esp(0)',
  'esp(1)',
  'esp(2)',
  'nsp(1)',
  'nsp(2)',
  'gwesp(0.05,fixed=TRUE)',
  'gwesp(0.25,fixed=TRUE)',
  'gwesp(0.5,fixed=TRUE)',
  'gwesp(0.75,fixed=TRUE)',
  'gwesp(2,fixed=TRUE)',
  'gwesp(1,fixed=TRUE)',
  'degcor'
)


#------------------------------------------------------------------------------
#NOTE - NOTHING BELOW THIS POINT SHOULD REQUIRE CHANGING UNLESS TINKERING WITH
#THE TECHNIQUE-----------------------------------------------------------------
#  (Should you do that? C'est vous qui voyez...)
#------------------------------------------------------------------------------


#LIBRARIES AND IMPORTS---------------------------------------------------------

#Load required libraries
library(sna)
library(ergm)
library(ergm.components)
library(ergm.multi)
library(parallel)           #Only needed if cores>1

#Load model fitting/evaluation functions
source("csh_cssc_1node_functions.R")


#LOGGING SETUP-----------------------------------------------------------------

# Setup logging to file
if(!dir.exists(logdir)) dir.create(logdir, recursive=TRUE)
if(!dir.exists(batch_log_dir)) dir.create(batch_log_dir, recursive=TRUE)

# Function to log both to console and file
log_message <- function(msg, file=logfile) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  cat(full_msg, "\n", file=file, append=TRUE)
}

# Function to redirect output to both console and log file
sink_to_log <- function(file=logfile) {
  # Create a connection that writes to both console and file
  con <- file(file, open="a")
  sink(con, append=TRUE, type="output")
  sink(con, append=TRUE, type="message")
}

# Function to stop redirecting output
stop_sink <- function() {
  sink(type="message")
  sink(type="output")
}

# Initialize log file
cat("", file=logfile)  # Clear the log file
log_message("=== RestFiber NHM Analysis Started ===")
log_message(paste("Base seed:", base.seed))
log_message(paste("Data directory:", datdir))
log_message(paste("Output directory:", outdir))
log_message(paste("Model file:", modfile))
log_message(paste("Checkpoint file:", checkpoint_file))
log_message(paste("Log file:", logfile))

#==============================================================================
#CHECKPOINT RECOVERY SYSTEM
#==============================================================================

# Check if both checkpoint and intermediate files exist - if so, skip all data processing
if(file.exists(checkpoint_file) && file.exists(intermediate_file)) {
  log_message("=== FULL CHECKPOINT RECOVERY MODE ===")
  log_message(paste("Found checkpoint file:", checkpoint_file))
  log_message(paste("Found intermediate file:", intermediate_file))
  log_message("Loading all previous progress and skipping data processing...")
  
  # Load checkpoint data
  load(checkpoint_file)
  
  # Load intermediate data
  load(intermediate_file)
  
  log_message("Successfully loaded checkpoint and intermediate data")
  log_message("Skipping data processing, baseline model fitting, and going straight to batch processing")
  
  # Set flags to skip everything except batch processing
  skip_data_processing <- TRUE
  skip_baseline_model <- TRUE
  skip_to_batch_processing <- TRUE
  
} else {

# Check if checkpoint file exists and load from there
if(file.exists(checkpoint_file)) {
  log_message("=== CHECKPOINT RECOVERY MODE ===")
  log_message(paste("Found checkpoint file:", checkpoint_file))
  log_message("Loading previous progress...")
  
  # Load checkpoint data
  load(checkpoint_file)
  
  # Check if mods contains valid models
  if(exists("mods")) {
    valid_mods_logical <- sapply(mods, function(x) !is.null(x) && is.list(x))
    valid_mods_count <- sum(valid_mods_logical)
    total_mods_count <- length(mods)
    log_message(paste("Found", valid_mods_count, "valid models out of", total_mods_count, "in checkpoint"))
    log_message(paste("Valid models (by index):", paste(which(valid_mods_logical), collapse = ", ")))
    if(valid_mods_count < total_mods_count) {
      log_message("WARNING: Checkpoint contains at least one invalid model - restarting batch processing")
      log_message("Clearing invalid mods variable and restarting batch processing...")
      # Clear the invalid mods variable but keep other valuable data
      mods <- NULL
      skip_to_model_selection <- FALSE
      skip_to_batch_processing <- FALSE
    } else {
      log_message("Checkpoint contains valid models - proceeding to model selection")
      skip_to_model_selection <- TRUE
    }
  } else {
    log_message("No mods variable in checkpoint - restarting batch processing")
    skip_to_model_selection <- FALSE
  }
  
} else if(file.exists(batch_checkpoint_file)) {
  log_message("=== BATCH CHECKPOINT RECOVERY MODE ===")
  log_message(paste("Found batch checkpoint file:", batch_checkpoint_file))
  log_message("Loading batch processing progress...")
  
  # Load batch checkpoint data
  load(batch_checkpoint_file)
  
  # Check if candidate list is identical
  if(exists("checkpoint_cand") && identical(cand, checkpoint_cand)) {
    log_message("Candidate list is identical - resuming from batch checkpoint")
    log_message(paste("Recovered progress: tested", term_counter, "terms"))
    log_message(paste("Recovered batch progress: completed", completed_batches, "batches"))
    log_message("Skipping data processing and resuming batch processing...")
    
    # Set flag to skip to batch processing
    skip_to_batch_processing <- TRUE
    skip_to_model_selection <- FALSE
  } else {
    log_message("Candidate list has changed - starting fresh analysis")
    log_message("Removing outdated batch checkpoint...")
    file.remove(batch_checkpoint_file)
    skip_to_batch_processing <- FALSE
    skip_to_model_selection <- FALSE
  }
  
} else {
  log_message("No checkpoint file found - starting fresh analysis")
  skip_to_model_selection <- FALSE
  skip_to_batch_processing <- FALSE
}

} # End of checkpoint recovery system

# Start redirecting output to log file
sink_to_log()

#==============================================================================
#MAIN ANALYSIS - CONDITIONAL EXECUTION BASED ON CHECKPOINT
#==============================================================================

if(!skip_data_processing && !skip_to_model_selection && !skip_to_batch_processing) {
  log_message("=== STARTING FRESH ANALYSIS ===")

#DATA PROCESSING---------------------------------------------------------------

#Load and process the restrained fiber data
log_message("Loading restrained fiber edgelist data...")

# Function to convert edgelist to network objects with proper node attributes
convert_edgelist_to_networks <- function(edgelist_file, system_type) {
  log_message(paste("Processing", system_type, "system from", basename(edgelist_file)))
  
  # Load the edgelist data
  load(edgelist_file)  # This loads 'edgelist' object
  
  networks <- list()
  
  for(i in 1:length(edgelist)) {
    # Get the edgelist matrix (3 columns: node1, node2, weight)
    el_matrix <- edgelist[[i]]
    
    # Determine total number of nodes
    max_node <- max(c(el_matrix[,1], el_matrix[,2]))
    
    # Create network from edgelist (using first 2 columns only)
    net <- network(el_matrix[,1:2], directed=FALSE, vertices=1:max_node)
    
    # Add node attributes based on system type
    if(system_type == "mix") {
      # Mix system: first 411 are CSSC, rest are CSH
      resnames <- c(rep("CSSC", 411), rep("CSH", max_node - 411))
    } else {
      # Pure system: all CSSC
      resnames <- rep("CSSC", max_node)
    }
    
    # Set vertex attributes
    net %v% "Index" <- 1:max_node
    net %v% "Resid" <- 1:max_node  # Using node index as residue ID
    net %v% "Resname" <- resnames
    
    networks[[i]] <- net
  }
  
  log_message(paste("Converted", length(networks), "networks for", system_type, "system"))
  return(networks)
}

# Process both system types
log_message("Converting mix system (CSSC + CSH)...")
mix_networks <- convert_edgelist_to_networks(mix_file, "mix")

log_message("Converting pure system (CSSC only)...")
pure_networks <- convert_edgelist_to_networks(pure_file, "pure")

# Create the nets.1n.cg structure expected by the rest of the script
# Index 1: Mix system (~50% CSSC based on 411/total nodes)  
# Index 2: Pure CSSC (100% CSSC)
nets.1n.cg <- list()
nets.1n.cg[[1]] <- mix_networks    # ~50% CSSC
nets.1n.cg[[2]] <- pure_networks   # 100% CSSC

log_message(paste("Created", length(nets.1n.cg), "compositional conditions"))
log_message(paste("Condition 1 (Mix):", length(nets.1n.cg[[1]]), "networks"))
log_message(paste("Condition 2 (Pure CSSC):", length(nets.1n.cg[[2]]), "networks"))

#DYNAMIC SAMPLE SIZE CALCULATION--------------------------------------------------
log_message("Calculating dynamic sample sizes based on available data...")

# Get network counts per condition
condition_sizes <- sapply(nets.1n.cg, length)
min_networks <- min(condition_sizes)

log_message(paste("Network counts per condition:", paste(condition_sizes, collapse=", ")))
log_message(paste("Minimum networks across conditions:", min_networks))

# Calculate sample sizes based on minimum available data and ratios
sspercond.train <- floor(min_networks * train_ratio)
sspercond.test <- floor(min_networks * test_ratio)
ss.gof <- floor(min_networks * gof_ratio)

# Ensure we don't exceed available data (add small buffer)
total_needed <- sspercond.train + sspercond.test + ss.gof
if(total_needed > min_networks) {
  # Adjust proportionally if needed
  scale_factor <- (min_networks - 1) / total_needed  # -1 for safety buffer
  sspercond.train <- floor(sspercond.train * scale_factor)
  sspercond.test <- floor(sspercond.test * scale_factor)
  ss.gof <- floor(ss.gof * scale_factor)
  log_message("Sample sizes adjusted to fit within available data")
}

log_message(paste("Final sample sizes per condition:"))
log_message(paste("  Training:", sspercond.train))
log_message(paste("  Test:", sspercond.test))
log_message(paste("  GOF:", ss.gof))
log_message(paste("  Total per condition:", sspercond.train + sspercond.test + ss.gof, "out of", min_networks, "available"))
log_message(paste("  Utilization:", round(100 * (sspercond.train + sspercond.test + ss.gof) / min_networks, 1), "%"))

#FAILSAFE CHECKS----------------------------------------------------------------
cat("=== FAILSAFE CHECKS ===\n")

# Check available node types across all conditions
all_node_types <- character(0)
for(i in 1:length(nets.1n.cg)) {
  for(j in 1:min(5, length(nets.1n.cg[[i]]))) {  # Check first 5 networks per condition
    net <- nets.1n.cg[[i]][[j]]
    node_types <- net %v% "Resname"
    all_node_types <- c(all_node_types, unique(node_types))
  }
}
unique_node_types <- unique(all_node_types)
cat("Available node types in data:", paste(unique_node_types, collapse=", "), "\n")

# Check edge types across conditions
cat("\nChecking edge type availability across conditions:\n")
edge_type_summary <- list()
for(i in 1:length(nets.1n.cg)) {
  cat("Condition", i, ":\n")
  
  # Sample a few networks to check edge patterns
  edge_counts <- list()
  for(j in 1:min(3, length(nets.1n.cg[[i]]))) {
    net <- nets.1n.cg[[i]][[j]]
    node_types <- net %v% "Resname"
    adj_mat <- as.matrix(net)
    
    # Count different edge types
    for(type1 in unique_node_types) {
      for(type2 in unique_node_types) {
        if(type1 <= type2) {  # Avoid double counting in undirected networks
          nodes1 <- which(node_types == type1)
          nodes2 <- which(node_types == type2)
          
          if(length(nodes1) > 0 && length(nodes2) > 0) {
            if(type1 == type2) {
              # Same type edges (homophily)
              if(length(nodes1) >= 2) {
                edge_count <- sum(adj_mat[nodes1, nodes1]) / 2
                edge_type <- paste0(type1, "-", type2)
              } else {
                edge_count <- 0
                edge_type <- paste0(type1, "-", type2)
              }
            } else {
              # Different type edges (heterophily)
              edge_count <- sum(adj_mat[nodes1, nodes2])
              edge_type <- paste0(type1, "-", type2)
            }
            
            if(is.null(edge_counts[[edge_type]])) {
              edge_counts[[edge_type]] <- edge_count
            } else {
              edge_counts[[edge_type]] <- edge_counts[[edge_type]] + edge_count
            }
          }
        }
      }
    }
  }
  
  # Print edge type availability for this condition
  for(edge_type in names(edge_counts)) {
    cat("  ", edge_type, "edges:", edge_counts[[edge_type]], "\n")
  }
  edge_type_summary[[i]] <- edge_counts
}

# Validate resnames configuration
cat("\nValidating resnames configuration:\n")
cat("Script expects resnames:", paste(resnames, collapse=", "), "\n")
cat("Data contains node types:", paste(unique_node_types, collapse=", "), "\n")

# Check if all expected resnames are present in data
missing_types <- setdiff(resnames, unique_node_types)
extra_types <- setdiff(unique_node_types, resnames)

if(length(missing_types) > 0) {
  cat("ERROR: Expected node types not found in data:", paste(missing_types, collapse=", "), "\n")
  cat("ABORTING ANALYSIS\n")
  stop("Resnames configuration mismatch - missing expected node types")
}

if(length(extra_types) > 0) {
  cat("WARNING: Additional node types found in data:", paste(extra_types, collapse=", "), "\n")
  cat("These will be ignored in analysis\n")
}

cat("SUCCESS: Resnames configuration is compatible with data\n")
cat("Proceeding with analysis using resnames:", paste(resnames, collapse=", "), "\n")

# Check if key edge types are available for nodematch terms
cat("\nChecking nodematch term viability:\n")
for(resname in resnames) {
  homophily_edges <- 0
  for(i in 1:length(edge_type_summary)) {
    edge_type <- paste0(resname, "-", resname)
    if(!is.null(edge_type_summary[[i]][[edge_type]])) {
      homophily_edges <- homophily_edges + edge_type_summary[[i]][[edge_type]]
    }
  }
  cat("  nodematch.Resname.", resname, " viable: ", homophily_edges, " total ", resname, "-", resname, " edges\n", sep="")
  if(homophily_edges == 0) {
    cat("    WARNING: No ", resname, "-", resname, " edges found - nodematch term may be unidentifiable\n", sep="")
  }
}

cat("=== END FAILSAFE CHECKS ===\n\n")

#Determine maximum possible number of neighbors; we use the max ever observed
#as a good proxy (given the very large number of draws, many involving compact
#phases)
#Add safety margin to avoid degree bound violations
maxdegs<-sapply(nets.1n.cg,function(z){
  dm<-sapply(z,function(w){
    d<-degree(w,gmode="graph")
    type<-w%v%"Resname"
    if(any(type=="CSH"))
      dm<-max(d[type=="CSH"])
    else
      dm<-0
    if(any(type=="CSSC"))
      dm<-c(dm,max(d[type=="CSSC"]))
    else
      dm<-c(dm,0)
    dm
  })
  apply(dm,1,max)
})
dmax<-apply(maxdegs,1,max)
#Add larger safety margin to prevent degree bound violations in restrained fiber data
dmax<-dmax + 10
cat("Maximum degrees calculated:", dmax, "\n")
cat("Degree bounds will be set to:", dmax, "\n")


#Merge everything into a master list - but keep it balanced; given the sample
#sizes, we draw the specified number of random observations from each condition
set.seed(base.seed)
trainsel<-lapply(1:length(nets.1n.cg),function(z){sample(1:length(nets.1n.cg[[z]]),sspercond.train)})
testsel<-lapply(1:length(nets.1n.cg),function(z){sample((1:length(nets.1n.cg[[z]]))[-trainsel[[z]]], sspercond.test)})
nets.train<-list()
for(i in 1:length(nets.1n.cg))
  nets.train<-c(nets.train,nets.1n.cg[[i]][trainsel[[i]]])
class(nets.train)<-"network.list"
nets.test<-list()
for(i in 1:length(nets.1n.cg))
  nets.test<-c(nets.test,nets.1n.cg[[i]][testsel[[i]]])
class(nets.test)<-"network.list"

} # End of data processing section

#MODEL SELECTION AND INFERENCE-------------------------------------------------

#Seek an optimal model by forward selection, using the held-out deviance
seed<-base.seed
options(ergm.loglik.warn_dyads=FALSE)  #This is annoying and irrelevant

# Initialize periodic saving counter
term_counter <- 0

if(!skip_baseline_model) {
  log_message("=== FITTING BASELINE MODEL ===")
  best.fit<-ergmPool(NULL, train=nets.train, test=nets.test, dmax=dmax, thin=thin, subphases=subphases, seed=seed, a1=a1)
  best.dev<-c(best.fit$deviance.test,best.fit$deviance.test.se)
  best.terms<-rep(FALSE,length(cand))
  hist.dev<-best.dev                   #Keep a record of what we've seen
  hist.terms<-best.terms
  maxsedev<-se.safety*best.dev[2]      #Set the threshold for auto-rejection due to low precision
  flag<-TRUE
  terminc<-best.terms
  isgwesp<-grepl("^gwesp",cand)        #Mark the curved terms - need to make mutually exclusive
  isgwdegree<-grepl("^gwdegree",cand)

  # Initial save of baseline model
  save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, flag, isgwesp, isgwdegree, maxsedev, file=intermediate_file)
  cat("Saved initial baseline model to", intermediate_file, "\n")
} else {
  log_message("=== SKIPPING BASELINE MODEL FITTING ===")
  log_message("Baseline model already exists in checkpoint/intermediate files")
  log_message("Initializing variables from checkpoint/intermediate data...")
  
  # Initialize variables that would normally be set during baseline model fitting
  # These should be loaded from checkpoint/intermediate files, but ensure they exist
  if(!exists("flag")) flag <- TRUE
  if(!exists("terminc")) terminc <- best.terms
  if(!exists("isgwesp")) isgwesp <- grepl("^gwesp",cand)
  if(!exists("isgwdegree")) isgwdegree <- grepl("^gwdegree",cand)
  if(!exists("maxsedev")) maxsedev <- se.safety*best.dev[2]
  
  log_message("Variables initialized successfully")
}
# Handle batch processing recovery
if(skip_to_batch_processing) {
  log_message("=== RESUMING BATCH PROCESSING ===")
  log_message("Skipped data processing and baseline model fitting - resuming from checkpoint")
  log_message("Proceeding directly to batch processing with existing data")
} else {
  log_message("=== STARTING BATCH PROCESSING ===")
}

while(flag&&any(!terminc)){
  flag<-FALSE
  if(cores>1){  #Parallel approach (be sure you have enough RAM)
    #Process candidates in batches to avoid resource exhaustion
    remaining_candidates <- which(!terminc)
    batch_size <- min(cores, length(remaining_candidates))
    
    # Initialize batch tracking variables
    if(!exists("completed_batches")) {
      completed_batches <- 0
    }
    if(!exists("mods")) {
      mods <- vector("list", length(cand))
    }
    
    cat("Processing", length(remaining_candidates), "remaining candidates in batches of", batch_size, "\n")
    cat("Starting from batch", completed_batches + 1, "\n")
    
    #Process candidates in batches
    for(batch_start in seq(1, length(remaining_candidates), by=batch_size)) {
      batch_end <- min(batch_start + batch_size - 1, length(remaining_candidates))
      batch_candidates <- remaining_candidates[batch_start:batch_end]
      current_batch <- ceiling(batch_start/batch_size)
      
      # Initialize batch log
      batch_log <- c()
      batch_log <- c(batch_log, paste("Processing batch", current_batch, "of", ceiling(length(remaining_candidates)/batch_size), 
                                     "- candidates", batch_start, "to", batch_end))
      
      stimecyc<-proc.time()[3]
      batch_mods<-mclapply(batch_candidates,function(i){
        #Set seed in a reproducible way
        locseed<-base.seed+i
        #Set up the candidate
        oterminc<-terminc
        terminc[i]<-TRUE
        if(isgwesp[i]){  #Only one gwesp at a time...
          terminc[oterminc&isgwesp]<-FALSE
        }else if(isgwdegree[i]){ #Only one gwdegree, too
          terminc[oterminc&isgwdegree]<-FALSE
        }
        # Store log message instead of printing
        work_msg <- paste("Working on",paste(cand[terminc],collapse="+"))
        stime<-proc.time()[3]
        #Fit the model with verbose=FALSE to avoid parallel output conflicts
        fit<-ergmPool(cand[terminc], train=nets.train, test=nets.test, dmax=dmax, thin=thin, seed=locseed, subphases=subphases, a1=a1, verbose=FALSE)
        finish_msg <- paste("\tFinished working on",paste(cand[terminc],collapse="+"),"- ET",(proc.time()[3]-stime)/60,"min")
        fit$terminc<-terminc
        fit$term_counter <- term_counter + 1  # Track which term this is
        fit$work_msg <- work_msg  # Store log messages in the result
        fit$finish_msg <- finish_msg
        #Return the result
        fit
      },mc.cores=cores, mc.preschedule=FALSE)
      
      #Store batch results in the main mods list and collect log messages
      batch_valid_count <- 0
      for(j in 1:length(batch_candidates)) {
        mods[[batch_candidates[j]]] <- batch_mods[[j]]
        # Collect log messages from the results
        if(!is.null(batch_mods[[j]]) && is.list(batch_mods[[j]])) {
          batch_valid_count <- batch_valid_count + 1
          if(!is.null(batch_mods[[j]]$work_msg)) {
            batch_log <- c(batch_log, batch_mods[[j]]$work_msg)
          }
          if(!is.null(batch_mods[[j]]$finish_msg)) {
            batch_log <- c(batch_log, batch_mods[[j]]$finish_msg)
          }
        } else {
          batch_log <- c(batch_log, paste("ERROR: Model", batch_candidates[j], "failed -", class(batch_mods[[j]])))
        }
      }
      
      # Check if batch had too many failures and retry if needed
      if(batch_valid_count < length(batch_candidates) * 0.5) {  # If less than 50% valid
        batch_log <- c(batch_log, paste("WARNING: Only", batch_valid_count, "out of", length(batch_candidates), "models valid in batch", current_batch))
        batch_log <- c(batch_log, "Retrying batch with reduced cores...")
        
        # Retry with fewer cores or single core
        retry_cores <- max(1, cores - 2)
        batch_log <- c(batch_log, paste("Retrying batch", current_batch, "with", retry_cores, "cores"))
        
        stimecyc_retry <- proc.time()[3]
        batch_mods_retry <- mclapply(batch_candidates, function(i) {
          #Set seed in a reproducible way
          locseed <- base.seed + i
          #Set up the candidate
          oterminc <- terminc
          terminc[i] <- TRUE
          if(isgwesp[i]) {
            terminc[oterminc & isgwesp] <- FALSE
          } else if(isgwdegree[i]) {
            terminc[oterminc & isgwdegree] <- FALSE
          }
          # Store log message instead of printing
          work_msg <- paste("RETRY: Working on", paste(cand[terminc], collapse="+"))
          stime <- proc.time()[3]
          #Fit the model with verbose=FALSE to avoid parallel output conflicts
          fit <- ergmPool(cand[terminc], train=nets.train, test=nets.test, dmax=dmax, thin=thin, seed=locseed, subphases=subphases, a1=a1, verbose=FALSE)
          finish_msg <- paste("\tRETRY: Finished working on", paste(cand[terminc], collapse="+"), "- ET", (proc.time()[3]-stime)/60, "min")
          fit$terminc <- terminc
          fit$term_counter <- term_counter + 1
          fit$work_msg <- work_msg
          fit$finish_msg <- finish_msg
          fit
        }, mc.cores=retry_cores, mc.preschedule=FALSE)
        
        # Update results with retry
        batch_valid_count_retry <- 0
        for(j in 1:length(batch_candidates)) {
          mods[[batch_candidates[j]]] <- batch_mods_retry[[j]]
          if(!is.null(batch_mods_retry[[j]]) && is.list(batch_mods_retry[[j]])) {
            batch_valid_count_retry <- batch_valid_count_retry + 1
            if(!is.null(batch_mods_retry[[j]]$work_msg)) {
              batch_log <- c(batch_log, batch_mods_retry[[j]]$work_msg)
            }
            if(!is.null(batch_mods_retry[[j]]$finish_msg)) {
              batch_log <- c(batch_log, batch_mods_retry[[j]]$finish_msg)
            }
          } else {
            batch_log <- c(batch_log, paste("RETRY ERROR: Model", batch_candidates[j], "still failed -", class(batch_mods_retry[[j]])))
          }
        }
        
        batch_log <- c(batch_log, paste("Retry completed:", batch_valid_count_retry, "out of", length(batch_candidates), "models valid"))
        batch_log <- c(batch_log, paste("Retry took", (proc.time()[3] - stimecyc_retry)/60, "min"))
      }
      
      completed_batches <- completed_batches + 1
      batch_log <- c(batch_log, paste("Finished batch", current_batch, "- took", (proc.time()[3]-stimecyc)/60, "min"))
      
      #Save batch checkpoint after each batch
      checkpoint_cand <- cand  # Save candidate list for comparison
      save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, mods, completed_batches, checkpoint_cand, flag, isgwesp, isgwdegree, maxsedev, file=batch_checkpoint_file)
      batch_log <- c(batch_log, paste("Saved batch checkpoint after batch", current_batch))
      
      # Write batch log to file
      batch_log_file <- paste0(batch_log_dir, "/batch_", current_batch, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
      writeLines(batch_log, batch_log_file)
      
      # Also log to main log file (single-threaded)
      log_message(paste("Batch", current_batch, "completed - log saved to", basename(batch_log_file)))
    }
    
    cat("Finished with all candidates in this cycle - finding best model\n")
    
    # DIAGNOSTIC: Check mods structure before processing
    cat("DIAGNOSTIC: mods has", length(mods), "elements\n")
    valid_mods_count <- 0
    for(j in 1:length(mods)) {
      if(is.null(mods[[j]])) {
        cat("WARNING: mods[[", j, "]] is NULL for term:", cand[j], "\n")
      } else if(!is.list(mods[[j]])) {
        cat("WARNING: mods[[", j, "]] is not a list for term:", cand[j], ", class:", class(mods[[j]]), "\n")
      } else {
        cat("OK: mods[[", j, "]] is a proper list for term:", cand[j], "\n")
        valid_mods_count <- valid_mods_count + 1
      }
    }
    
    # Check if we have enough valid models to proceed
    if(valid_mods_count == 0) {
      cat("ERROR: No valid models found in this cycle!\n")
      cat("This suggests a systematic issue with the model fitting process.\n")
      cat("Consider reducing cores, checking data quality, or adjusting model parameters.\n")
      # Continue anyway to avoid infinite loop, but flag will remain FALSE
    } else {
      cat("Found", valid_mods_count, "valid models out of", length(mods), "total\n")
    }
    #See if any of these are better than what we have, keeping the best
    # CHECKPOINT SAVE before critical model selection loop
    checkpoint_before_selection <- paste0(outdir, "/csh_cssc_1node_fit_RestFiber_before_selection.Rdata")
    save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, mods, flag, isgwesp, isgwdegree, maxsedev, file=checkpoint_before_selection)
    cat("CHECKPOINT: Saved state before model selection to", checkpoint_before_selection, "\n")
    
    # FAILSAFE model selection loop with detailed error checking
    # Process ALL candidates that were tested in this cycle (not just remaining ones)
    tested_candidates <- which(!is.null(mods) & sapply(mods, function(x) !is.null(x)))
    cat("Evaluating", length(tested_candidates), "tested candidates for model selection\n")
    
    for(i in tested_candidates){
      term_counter <- term_counter + 1  # Increment counter for each term tested
      
      # FAILSAFE: Check model structure before accessing
      cat("Checking model", i, "for term:", cand[i], "\n")
      
      # Check if mods[[i]] is a proper list
      if(!is.list(mods[[i]])) {
        cat("ERROR: mods[[", i, "]] is not a list, it's a", class(mods[[i]]), "\n")
        cat("Term:", cand[i], "\n")
        cat("Value:", mods[[i]], "\n")
        next  # Skip this model and continue
      }
      
      # Check if required fields exist
      if(is.null(mods[[i]]$deviance.test.se)) {
        cat("ERROR: mods[[", i, "]]$deviance.test.se is NULL for term:", cand[i], "\n")
        cat("Available fields:", names(mods[[i]]), "\n")
        next  # Skip this model and continue
      }
      
      if(is.null(mods[[i]]$deviance.test)) {
        cat("ERROR: mods[[", i, "]]$deviance.test is NULL for term:", cand[i], "\n")
        cat("Available fields:", names(mods[[i]]), "\n")
        next  # Skip this model and continue
      }
      
      # Now safely access the fields
      tryCatch({
        if((mods[[i]]$deviance.test.se<maxsedev) && ((mods[[i]]$deviance.test-best.dev[1])/sqrt(mods[[i]]$deviance.test.se^2+best.dev[2]^2) < -se.accept)){
          cat("\tAdding term",cand[i],"was an improvement!\n")
          best.fit<-mods[[i]]
          best.terms<-mods[[i]]$terminc
          best.dev<-c(mods[[i]]$deviance.test,mods[[i]]$deviance.test.se)
          flag<-TRUE
        }
      }, error = function(e) {
        cat("ERROR in model comparison for i =", i, ", term =", cand[i], "\n")
        cat("Error message:", e$message, "\n")
        cat("mods[[i]] structure:\n")
        str(mods[[i]])
      })
    }
    
    # Periodic saving after testing terms (during model selection)
    if(term_counter %% save_interval == 0) {
      cat("Periodic save at term", term_counter, "\n")
      save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, flag, isgwesp, isgwdegree, maxsedev, file=intermediate_file)
      cat("Saved intermediate results to", intermediate_file, "\n")
    }
  }else{  #Non-parallel approach (e.g., for Windows)
    stimecyc<-proc.time()[3]
    # Initialize cycle log for non-parallel processing
    cycle_log <- c()
    cycle_log <- c(cycle_log, paste("Processing", length(which(!terminc)), "candidates sequentially"))
    
    for(i in which(!terminc)){
      #Set seed in a reproducible way
      locseed<-base.seed+i
      #Set up the candidate
      oterminc<-terminc
      terminc[i]<-TRUE
      if(isgwesp[i]){  #Only one gwesp at a time...
        terminc[oterminc&isgwesp]<-FALSE
      }else if(isgwdegree[i]){ #Only one gwdegree, too
        terminc[oterminc&isgwdegree]<-FALSE
      }
      
      work_msg <- paste("Working on",paste(cand[terminc],collapse="+"))
      cycle_log <- c(cycle_log, work_msg)
      cycle_log <- c(cycle_log, paste("Best deviance",best.dev[1],"+/-",best.dev[2]))
      
      stime<-proc.time()[3]
      #Fit the model with verbose=FALSE to avoid output conflicts
      fit<-ergmPool(cand[terminc], train=nets.train, test=nets.test, dmax=dmax, thin=thin, seed=locseed, subphases=subphases, a1=a1, verbose=FALSE)
      finish_msg <- paste("\tFinished working on",paste(cand[terminc],collapse="+"),"- ET",(proc.time()[3]-stime)/60,"min")
      cycle_log <- c(cycle_log, finish_msg)
      
      term_counter <- term_counter + 1  # Increment counter for each term tested
      
      if((fit$deviance.test.se<maxsedev) && ((fit$deviance.test-best.dev[1])/sqrt(fit$deviance.test.se^2+best.dev[2]^2) < -se.accept)){
        improvement_msg <- paste("\tAdding term",cand[i],"was an improvement!")
        cycle_log <- c(cycle_log, improvement_msg)
        best.fit<-fit
        best.terms<-terminc
        best.dev<-c(fit$deviance.test,fit$deviance.test.se)
        flag<-TRUE
      }
      
      # Periodic saving after testing each term
      if(term_counter %% save_interval == 0) {
        save_msg <- paste("Periodic save at term", term_counter)
        cycle_log <- c(cycle_log, save_msg)
        save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, flag, isgwesp, isgwdegree, maxsedev, file=intermediate_file)
        cycle_log <- c(cycle_log, paste("Saved intermediate results to", basename(intermediate_file)))
      }
      
      terminc<-oterminc
    }
    
    cycle_log <- c(cycle_log, "Finished with all candidates in this cycle - finding best model")
    cycle_log <- c(cycle_log, paste("\tCycle took",(proc.time()[3]-stimecyc)/60,"min"))
    
    # Write cycle log to file
    cycle_log_file <- paste0(batch_log_dir, "/cycle_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
    writeLines(cycle_log, cycle_log_file)
    log_message(paste("Cycle completed - log saved to", basename(cycle_log_file)))
  }
  terminc<-best.terms                        #Update inclusions
  if(flag){
    hist.terms<-rbind(hist.terms,best.terms)   #Update the history
    hist.dev<-rbind(hist.dev,best.dev)
  }
} # End of batch processing section

#==============================================================================
#FAILSAFE CHECKPOINT SAVE AFTER COMPLETING ALL CANDIDATE TESTING
#==============================================================================

log_message("=== CREATING FAILSAFE CHECKPOINT ===")
log_message("All candidate terms have been tested - saving checkpoint before model selection")

# Check if we have valid models before saving checkpoint
if(exists("mods")) {
  valid_mods_count <- sum(sapply(mods, function(x) !is.null(x) && is.list(x)))
  log_message(paste("Found", valid_mods_count, "valid models to save in checkpoint"))
  
  if(valid_mods_count == 0) {
    log_message("WARNING: No valid models found in this cycle")
    log_message("This suggests a systematic failure in batch processing")
    log_message("Restarting batch processing with reduced cores...")
    
    # Clear the invalid mods and restart batch processing
    mods <- NULL
    

    
    # The while loop will continue and retry batch processing
  } else {
    # Clean up batch checkpoint file since we're done with batch processing
    if(file.exists(batch_checkpoint_file)) {
      file.remove(batch_checkpoint_file)
      log_message("Removed batch checkpoint file - batch processing complete")
    }
    
    # Save all critical variables needed for model selection
    save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, 
         terminc, term_counter, mods, nets.train, nets.test, ss.gof, flag, isgwesp, isgwdegree, maxsedev,
         file=checkpoint_file)
    
    log_message(paste("Checkpoint saved to:", checkpoint_file))
    log_message("If analysis fails during model selection, restart will resume from here")
  }
} else {
  log_message("ERROR: No mods variable found - cannot create checkpoint")
}

} else {
  log_message("=== RESUMING FROM CHECKPOINT ===")
  log_message("Skipped data processing and candidate testing")
}

#==============================================================================
#MODEL SELECTION PHASE (ALWAYS EXECUTED)
#==============================================================================

log_message("=== STARTING MODEL SELECTION PHASE ===")

# Enhanced error handling for model selection
tryCatch({
  
  # Check if we have the mods variable (from fresh run or checkpoint)
  if(!exists("mods")) {
    log_message("ERROR: mods variable not found - cannot proceed with model selection")
    log_message("This suggests the candidate testing phase did not complete properly")
    stop("Missing mods variable required for model selection")
  }
  
  # Validate mods structure before proceeding
  log_message(paste("Found", length(mods), "candidate models to evaluate"))
  
  # Check each model in mods for proper structure
  valid_mods <- sapply(mods, function(mod) {
    is.list(mod) && !is.null(mod$deviance.test) && !is.null(mod$deviance.test.se)
  })
  
  if(!all(valid_mods)) {
    log_message(paste("WARNING:", sum(!valid_mods), "models have invalid structure"))
    log_message("Filtering out invalid models before selection...")
    mods <- mods[valid_mods]
  }
  
  log_message(paste("Proceeding with", length(mods), "valid models"))

#Save the model info (w/out adequacy checks, lest we get interrupted)
save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, term_counter, file=modfile)
cat("Final model saved after testing", term_counter, "terms total\n")


#ADEQUACY CHECKS---------------------------------------------------------------

#Perform some basic adequacy checks
set.seed(seed)
best.gof<-gof(best.fit, dat=nets.test, nets=ss.gof, reps=500, verbose=TRUE)

#Save the adequacy check results (we add them to the above)
save(seed, cand, maxdegs, dmax, best.fit, best.dev, best.terms, hist.dev, hist.terms, terminc, best.gof, term_counter, file=modfile)
cat("Final analysis complete! Tested", term_counter, "terms total\n")
cat("RestFiber analysis complete!\n")

}, error = function(e) {
  log_message("=== ERROR DURING MODEL SELECTION ===")
  log_message(paste("Error message:", e$message))
  log_message("Checkpoint file preserved for recovery")
  log_message(paste("To recover, restart the script - it will resume from:", checkpoint_file))
  stop(paste("Model selection failed:", e$message))
})

# Stop logging and close connections
stop_sink()

#==============================================================================
#CLEANUP CHECKPOINT FILE ON SUCCESSFUL COMPLETION
#==============================================================================

# Remove checkpoint file since analysis completed successfully
if(file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
  log_message("Checkpoint file removed - analysis completed successfully")
}

log_message("=== RestFiber NHM Analysis Completed Successfully ===")
log_message(paste("Total terms tested:", term_counter))
log_message(paste("Final model saved to:", modfile))
log_message(paste("Log file saved to:", logfile))
