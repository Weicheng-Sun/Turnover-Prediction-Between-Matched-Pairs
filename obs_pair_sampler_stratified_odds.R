##############################################################################
## Distance-Stratified Observation Pair Sampler
##
## target_nMatch:
##   Maximum final number of pairs retained per distance bin.
##
## The mismatch proportion is calculated separately for each distance bin
## using all unique pairs collected after distance filtering.
##############################################################################

library(geosphere)
library(data.table)


get_distance_bins <- function() {
  data.frame(
    bin   = 1:6,
    lower = c(0.5, 8, 16, 32, 64, 128),
    upper = c(8, 16, 32, 64, 128, 256),
    label = c("1-8", "8-16", "16-32", "32-64", "64-128", "128-256")
  )
}


assign_distance_bin <- function(
    dist_km,
    distance_bins = get_distance_bins()
) {
  idx <- rep(NA_integer_, length(dist_km))
  
  for (i in seq_len(nrow(distance_bins))) {
    idx[
      dist_km >= distance_bins$lower[i] &
        dist_km < distance_bins$upper[i]
    ] <- distance_bins$bin[i]
  }
  
  idx
}


calc_distance_km <- function(
    lon1,
    lat1,
    lon2,
    lat2
) {
  geosphere::distHaversine(
    cbind(lon1, lat1),
    cbind(lon2, lat2)
  ) / 1000
}


##############################################################################
## Construct an order-independent ID for each observation pair
##############################################################################

make_pair_id <- function(pair_data) {
  
  n <- nrow(pair_data)
  
  if (n == 0L) {
    return(character(0))
  }
  
  if ("observation.numbers" %in% names(pair_data)) {
    
    obs <- pair_data[["observation.numbers"]]
    
    # Case 1: two-column matrix or data.frame
    if (is.matrix(obs) || is.data.frame(obs)) {
      if (ncol(obs) == 2L) {
        id1 <- as.character(obs[, 1])
        id2 <- as.character(obs[, 2])
        
        return(paste(
          pmin(id1, id2),
          pmax(id1, id2),
          sep = "~"
        ))
      }
    }
    
    # Case 2: list-column in which each row contains two IDs
    if (is.list(obs) && length(obs) == n && all(lengths(obs) == 2L)) {
      id1 <- vapply(
        obs,
        function(x) as.character(x[[1]]),
        character(1)
      )
      
      id2 <- vapply(
        obs,
        function(x) as.character(x[[2]]),
        character(1)
      )
      
      return(paste(
        pmin(id1, id2),
        pmax(id1, id2),
        sep = "~"
      ))
    }
    
    # Case 3: one already-constructed pair ID per row
    if (is.atomic(obs) && length(obs) == n) {
      return(as.character(obs))
    }
    
    # Case 4: flattened vector containing two IDs per row
    if (is.atomic(obs) && length(obs) == 2L * n) {
      obs_matrix <- matrix(
        obs,
        nrow = n,
        ncol = 2
      )
      
      id1 <- as.character(obs_matrix[, 1])
      id2 <- as.character(obs_matrix[, 2])
      
      return(paste(
        pmin(id1, id2),
        pmax(id1, id2),
        sep = "~"
      ))
    }
  }
  
  # Try commonly used pairs of observation-ID columns
  possible_id_columns <- list(
    c("observation.number1", "observation.number2"),
    c("observation.number.1", "observation.number.2"),
    c("observation_number_1", "observation_number_2"),
    c("observation_1", "observation_2"),
    c("Observation1", "Observation2"),
    c("Obs1", "Obs2"),
    c("ID1", "ID2")
  )
  
  for (id_columns in possible_id_columns) {
    if (all(id_columns %in% names(pair_data))) {
      id1 <- as.character(pair_data[[id_columns[1]]])
      id2 <- as.character(pair_data[[id_columns[2]]])
      
      return(paste(
        pmin(id1, id2),
        pmax(id1, id2),
        sep = "~"
      ))
    }
  }
  
  # Final fallback: construct the pair ID from the two endpoints
  coordinate_columns <- c(
    "Lon1", "Lat1",
    "Lon2", "Lat2"
  )
  
  if (all(coordinate_columns %in% names(pair_data))) {
    
    point1 <- paste(
      format(pair_data$Lon1, digits = 15, scientific = FALSE, trim = TRUE),
      format(pair_data$Lat1, digits = 15, scientific = FALSE, trim = TRUE),
      sep = ","
    )
    
    point2 <- paste(
      format(pair_data$Lon2, digits = 15, scientific = FALSE, trim = TRUE),
      format(pair_data$Lat2, digits = 15, scientific = FALSE, trim = TRUE),
      sep = ","
    )
    
    return(paste(
      pmin(point1, point2),
      pmax(point1, point2),
      sep = "~"
    ))
  }
  
  stop(
    paste0(
      "Unable to construct pair IDs. The sampler output must contain either ",
      "observation IDs or Lon1, Lat1, Lon2 and Lat2."
    )
  )
}


##############################################################################
## Main sampling function
##############################################################################

obsPairSampler.distanceStratified_2 <- function(
    frog.auGrid,
    m1,
    nMatch,
    target_nMatch = NULL,
    richness = TRUE,
    speciesThreshold = 500,
    coresToUse = parallel::detectCores() - 1,
    max_iterations = 20,
    distance_bins = get_distance_bins(),
    seed = NULL
) {
  
  if (is.null(target_nMatch)) {
    stop("Please specify target_nMatch.")
  }
  
  if (
    length(target_nMatch) != 1L ||
    !is.finite(target_nMatch) ||
    target_nMatch <= 0
  ) {
    stop("target_nMatch must be one positive number.")
  }
  
  target_nMatch <- as.integer(target_nMatch)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  cat(sprintf(
    paste0(
      "Target per distance bin: up to %d total pairs\n",
      "The match:mismatch ratio will be calculated separately ",
      "for each distance bin.\n"
    ),
    target_nMatch
  ))
  
  all_pairs <- NULL
  all_pair_ids <- character(0)
  
  ##########################################################################
  ## Sampling loop
  ##########################################################################
  
  for (iter in seq_len(max_iterations)) {
    
    cat(sprintf(
      "\nIteration %d/%d\n",
      iter,
      max_iterations
    ))
    
    new_pairs <- obsPairSampler.bigData.RECA(
      frog.auGrid      = frog.auGrid,
      nMatch           = nMatch,
      m1               = m1,
      richness         = richness,
      speciesThreshold = speciesThreshold,
      coresToUse       = coresToUse
    )
    
    if (is.null(new_pairs) || nrow(new_pairs) == 0L) {
      cat("No pairs returned in this iteration.\n")
      next
    }
    
    ########################################################################
    ## Calculate distance and retain pairs between 0.5 and 256 km
    ########################################################################
    
    dist_km <- calc_distance_km(
      new_pairs$Lon1,
      new_pairs$Lat1,
      new_pairs$Lon2,
      new_pairs$Lat2
    )
    
    bin_idx <- assign_distance_bin(
      dist_km,
      distance_bins
    )
    
    keep <- !is.na(bin_idx)
    
    new_pairs <- new_pairs[
      keep,
      ,
      drop = FALSE
    ]
    
    if (nrow(new_pairs) == 0L) {
      cat("No pairs within the specified distance range.\n")
      next
    }
    
    new_pairs$Distance_km  <- dist_km[keep]
    new_pairs$Distance_bin <- bin_idx[keep]
    
    ########################################################################
    ## Remove duplicate pairs within the current iteration
    ########################################################################
    
    new_pair_ids <- make_pair_id(new_pairs)
    
    keep_unique <- !duplicated(new_pair_ids)
    
    new_pairs <- new_pairs[
      keep_unique,
      ,
      drop = FALSE
    ]
    
    new_pair_ids <- new_pair_ids[keep_unique]
    
    ########################################################################
    ## Remove pairs already collected in previous iterations
    ########################################################################
    
    if (length(all_pair_ids) > 0L) {
      keep_new <- !(new_pair_ids %in% all_pair_ids)
      
      new_pairs <- new_pairs[
        keep_new,
        ,
        drop = FALSE
      ]
      
      new_pair_ids <- new_pair_ids[keep_new]
    }
    
    if (nrow(new_pairs) == 0L) {
      cat("No new unique pairs in this iteration.\n")
      next
    }
    
    ########################################################################
    ## Accumulate unique pairs
    ########################################################################
    
    if (is.null(all_pairs)) {
      all_pairs <- new_pairs
    } else {
      all_pairs <- rbind(
        all_pairs,
        new_pairs
      )
    }
    
    all_pair_ids <- c(
      all_pair_ids,
      new_pair_ids
    )
    
    cat(sprintf(
      "Added %d unique pairs; accumulated %d unique pairs.\n",
      nrow(new_pairs),
      nrow(all_pairs)
    ))
    
    gc()
  }
  
  if (is.null(all_pairs) || nrow(all_pairs) == 0L) {
    stop("No pairs were collected.")
  }
  
  ##########################################################################
  ## Final duplicate check
  ##########################################################################
  
  final_pair_ids <- make_pair_id(all_pairs)
  
  keep_unique <- !duplicated(final_pair_ids)
  
  all_pairs <- all_pairs[
    keep_unique,
    ,
    drop = FALSE
  ]
  
  cat(sprintf(
    "\nUnique pairs after all iterations: %d\n",
    nrow(all_pairs)
  ))
  
  ##########################################################################
  ## Check Match variable
  ##########################################################################
  
  if (!"Match" %in% names(all_pairs)) {
    stop("The sampler output does not contain a Match column.")
  }
  
  if (anyNA(all_pairs$Match)) {
    stop("The Match column contains missing values.")
  }
  
  if (!all(all_pairs$Match %in% c(0, 1))) {
    stop("The Match column must contain only 0 and 1.")
  }
  
  ##########################################################################
  ## Calculate mismatch proportion and sample within each distance bin
  ##########################################################################
  
  final_list <- vector(
    "list",
    nrow(distance_bins)
  )
  
  sampling_df <- data.frame(
    bin                       = distance_bins$label,
    before_total              = integer(nrow(distance_bins)),
    before_match              = integer(nrow(distance_bins)),
    before_mismatch           = integer(nrow(distance_bins)),
    observed_mismatch_prop    = rep(NA_real_, nrow(distance_bins)),
    target_total              = integer(nrow(distance_bins)),
    target_match              = integer(nrow(distance_bins)),
    target_mismatch           = integer(nrow(distance_bins)),
    final_total               = integer(nrow(distance_bins)),
    final_match               = integer(nrow(distance_bins)),
    final_mismatch            = integer(nrow(distance_bins)),
    final_mismatch_prop       = rep(NA_real_, nrow(distance_bins)),
    stringsAsFactors          = FALSE
  )
  
  for (b in seq_len(nrow(distance_bins))) {
    
    bin_pairs <- all_pairs[
      all_pairs$Distance_bin == distance_bins$bin[b],
      ,
      drop = FALSE
    ]
    
    n_total <- nrow(bin_pairs)
    
    if (n_total == 0L) {
      next
    }
    
    match_pairs <- bin_pairs[
      bin_pairs$Match == 0,
      ,
      drop = FALSE
    ]
    
    mismatch_pairs <- bin_pairs[
      bin_pairs$Match == 1,
      ,
      drop = FALSE
    ]
    
    n_match    <- nrow(match_pairs)
    n_mismatch <- nrow(mismatch_pairs)
    
    mismatch_prop <- n_mismatch / n_total
    
    sampling_df$before_total[b]           <- n_total
    sampling_df$before_match[b]           <- n_match
    sampling_df$before_mismatch[b]        <- n_mismatch
    sampling_df$observed_mismatch_prop[b] <- mismatch_prop
    
    ########################################################################
    ## If fewer than target_nMatch pairs are available, retain all pairs
    ########################################################################
    
    target_total <- min(
      n_total,
      target_nMatch
    )
    
    ########################################################################
    ## Calculate target counts from the observed mismatch proportion
    ########################################################################
    
    if (target_total == n_total) {
      
      target_match    <- n_match
      target_mismatch <- n_mismatch
      
    } else {
      
      target_mismatch <- as.integer(round(
        target_total * mismatch_prop
      ))
      
      target_match <- target_total - target_mismatch
      
      # Safeguards for very small groups
      if (target_mismatch > n_mismatch) {
        target_mismatch <- n_mismatch
        target_match <- target_total - target_mismatch
      }
      
      if (target_match > n_match) {
        target_match <- n_match
        target_mismatch <- target_total - target_match
      }
    }
    
    sampling_df$target_total[b]    <- target_match + target_mismatch
    sampling_df$target_match[b]    <- target_match
    sampling_df$target_mismatch[b] <- target_mismatch
    
    ########################################################################
    ## Randomly sample Match and Mismatch separately
    ########################################################################
    
    if (target_match < n_match) {
      match_pairs <- match_pairs[
        sample.int(
          n_match,
          target_match,
          replace = FALSE
        ),
        ,
        drop = FALSE
      ]
    }
    
    if (target_mismatch < n_mismatch) {
      mismatch_pairs <- mismatch_pairs[
        sample.int(
          n_mismatch,
          target_mismatch,
          replace = FALSE
        ),
        ,
        drop = FALSE
      ]
    }
    
    selected_pairs <- rbind(
      match_pairs,
      mismatch_pairs
    )
    
    if (nrow(selected_pairs) > 1L) {
      selected_pairs <- selected_pairs[
        sample.int(nrow(selected_pairs)),
        ,
        drop = FALSE
      ]
    }
    
    final_list[[b]] <- selected_pairs
    
    ########################################################################
    ## Save final counts
    ########################################################################
    
    final_match <- sum(
      selected_pairs$Match == 0
    )
    
    final_mismatch <- sum(
      selected_pairs$Match == 1
    )
    
    final_total <- nrow(selected_pairs)
    
    sampling_df$final_total[b]    <- final_total
    sampling_df$final_match[b]    <- final_match
    sampling_df$final_mismatch[b] <- final_mismatch
    
    if (final_total > 0L) {
      sampling_df$final_mismatch_prop[b] <-
        final_mismatch / final_total
    }
  }
  
  ##########################################################################
  ## Combine final samples
  ##########################################################################
  
  keep_list <- vapply(
    final_list,
    function(x) {
      !is.null(x) && nrow(x) > 0L
    },
    logical(1)
  )
  
  final_list <- final_list[keep_list]
  
  if (length(final_list) == 0L) {
    stop("No pairs were retained in the final sample.")
  }
  
  final_pairs <- do.call(
    rbind,
    final_list
  )
  
  rownames(final_pairs) <- NULL
  
  ##########################################################################
  ## Complete diagnostic table
  ##########################################################################
  
  sampling_df$observed_mismatch_pct <-
    100 * sampling_df$observed_mismatch_prop
  
  sampling_df$final_mismatch_pct <-
    100 * sampling_df$final_mismatch_prop
  
  sampling_df$completion_pct <-
    100 * sampling_df$final_total / target_nMatch
  
  bin_counts <- sampling_df$final_total
  names(bin_counts) <- sampling_df$bin
  
  cat("\nSampling summary:\n")
  
  print(
    sampling_df[
      ,
      c(
        "bin",
        "before_match",
        "before_mismatch",
        "observed_mismatch_pct",
        "target_match",
        "target_mismatch",
        "final_match",
        "final_mismatch",
        "final_mismatch_pct",
        "completion_pct"
      )
    ],
    row.names = FALSE,
    digits = 4
  )
  
  cat(sprintf(
    "\nTotal final pairs: %d\n",
    nrow(final_pairs)
  ))
  
  cat(sprintf(
    "Total final matches: %d\n",
    sum(final_pairs$Match == 0)
  ))
  
  cat(sprintf(
    "Total final mismatches: %d\n",
    sum(final_pairs$Match == 1)
  ))
  
  ##########################################################################
  ## Return results
  ##########################################################################
  
  result <- list(
    pairs                = final_pairs,
    all_unique_pairs     = all_pairs,
    bin_counts           = bin_counts,
    sampling_df          = sampling_df,
    balance_df           = sampling_df,
    total_matches        = sum(final_pairs$Match == 0),
    total_mismatches     = sum(final_pairs$Match == 1),
    target_pairs_per_bin = target_nMatch,
    iterations           = max_iterations,
    distance_bins        = distance_bins
  )
  
  class(result) <- "obsPairSampler.distanceStratified"
  
  return(result)
}