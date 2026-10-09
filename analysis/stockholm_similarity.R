#!/usr/bin/env Rscript
# stockholm_pairwise_similarity.R
#
# Parses a set of Stockholm multiple sequence alignment files and computes,
# for every pair of files, two similarity measures:
#
#   1. n_shared_matrix   - number of sequence identifiers found in both files
#   2. pct_overlap_matrix - mean percentage coordinate overlap of the hits
#                           shared between the two files
#
# Assumes sequence IDs follow the common Infernal/HMMER convention:
#   <accession>/<start>-<end>
# e.g. "NC_000913.3/145032-145091". IDs without this pattern are still
# counted for "shared sequences" but contribute no coordinate information.
#
# Dependencies: base R only. (An earlier version of this script used
# Biostrings::readRNAMultipleAlignment() for parsing, but its Stockholm
# reader is strict about block/row ordering and fails on real Infernal
# output - e.g. files with interleaved #=GR PP lines - with "alignment
# rows out of order". Since we only need the sequence ID strings, not the
# parsed alignment matrix, plain line-scanning is both simpler and more
# robust here.)
#
# Usage:
#   Rscript stockholm_pairwise_similarity.R <directory_of_stockholm_files> [output_prefix]
#
# Outputs (written to the working directory, or as CSVs if output_prefix given):
#   <prefix>_n_shared.csv
#   <prefix>_pct_overlap.csv

## ---- 1. Parsing ------------------------------------------------------

#' Extract the unique sequence identifiers from a Stockholm file.
#'
#' Stockholm alignments can wrap sequences across multiple blocks, and
#' annotation lines (#=GF, #=GS, #=GC, #=GR) can appear interleaved with
#' sequence lines in any order - this only looks at each line in isolation,
#' so block/row ordering doesn't matter.
#'
#' @param file path to a Stockholm (.sto/.stk) file
#' @return character vector of unique sequence IDs (first whitespace-delimited
#'   field of each sequence line)
read_stockholm_ids <- function(file) {
  lines <- readLines(file, warn = FALSE)

  # Drop blank lines, alignment terminator, and all annotation lines
  # (#=GF, #=GS, #=GC, #=GR, and the "# STOCKHOLM 1.0" header).
  is_annotation <- grepl("^\\s*(#|//\\s*$)", lines)
  is_blank      <- grepl("^\\s*$", lines)
  seq_lines     <- lines[!is_annotation & !is_blank]

  if (length(seq_lines) == 0) {
    return(character(0))
  }

  ids <- vapply(strsplit(trimws(seq_lines), "\\s+"), `[`, character(1), 1)
  unique(ids)
}

#' Parse sequence IDs into base accession + absolute coordinate range.
#'
#' Handles two ID conventions and resolves both down to a common
#' (accession, start, end) so hits from either naming scheme can be
#' compared against one another:
#'
#'   1. Simple:   <accession>/<start>-<end>
#'                e.g. "CP002116.1/581152-581028"
#'
#'   2. Extended: <cds_id>|<region_type>|<accession>:<region_start>-<region_end>(<strand>)/<rel_start>-<rel_end>
#'                e.g. "cds-AGW40642.1|utr5|CP006734.1:411858-413747(-)/1796-1890"
#'                Here <rel_start>-<rel_end> are coordinates *within* the
#'                extracted flanking region, not the genome, so they're
#'                converted to absolute genomic coordinates as
#'                  region_start + rel_pos - 1
#'
#' @param ids character vector of Stockholm sequence IDs
#' @return data.frame with columns id, base_id, start, end, strand.
#'   start/end are NA if the ID matches neither pattern.
parse_coordinates <- function(ids) {
  extended_pattern <- "^(.+)\\|(.+)\\|([^:]+):([0-9]+)-([0-9]+)\\(([+-])\\)/([0-9]+)-([0-9]+)$"
  simple_pattern    <- "^(.*)/([0-9]+)-([0-9]+)$"

  is_extended <- grepl(extended_pattern, ids)
  is_simple   <- !is_extended & grepl(simple_pattern, ids)

  base_id <- ids
  start   <- rep(NA_real_, length(ids))
  end     <- rep(NA_real_, length(ids))
  strand  <- rep(NA_character_, length(ids))

  # --- Extended CDS/UTR-relative IDs -----------------------------------
  # ASSUMPTION: rel_start/rel_end are numbered along the genome's forward
  # direction (region NOT reverse-complemented before numbering), even for
  # "(-)" strand entries - matching "region_start + rel_pos - 1". If your
  # pipeline instead numbers relative coordinates 5'->3' along the minus
  # strand, swap the two marked lines below for
  #   region_end[is_extended] - rel_start - 1   and
  #   region_end[is_extended] - rel_end   - 1
  if (any(is_extended)) {
    region_acc   <- sub(extended_pattern, "\\3", ids[is_extended])
    region_start <- as.numeric(sub(extended_pattern, "\\4", ids[is_extended]))
    rel_start    <- as.numeric(sub(extended_pattern, "\\7", ids[is_extended]))
    rel_end      <- as.numeric(sub(extended_pattern, "\\8", ids[is_extended]))

    base_id[is_extended] <- region_acc
    start[is_extended]   <- region_start + rel_start - 1  # <- flip here if minus-strand-numbered
    end[is_extended]     <- region_start + rel_end   - 1  # <- flip here if minus-strand-numbered
    strand[is_extended]  <- sub(extended_pattern, "\\6", ids[is_extended])
  }

  # --- Simple <accession>/<start>-<end> IDs ----------------------------
  if (any(is_simple)) {
    base_id[is_simple] <- sub(simple_pattern, "\\1", ids[is_simple])
    start[is_simple]   <- as.numeric(sub(simple_pattern, "\\2", ids[is_simple]))
    end[is_simple]     <- as.numeric(sub(simple_pattern, "\\3", ids[is_simple]))
  }

  # normalise so start <= end (reverse-strand hits can be reported as
  # start > end)
  lo <- pmin(start, end)
  hi <- pmax(start, end)

  data.frame(id = ids, base_id = base_id, start = lo, end = hi,
             strand = strand, stringsAsFactors = FALSE)
}

#' Load one Stockholm file into a parsed hit table.
load_stockholm_hits <- function(file) {
  parse_coordinates(read_stockholm_ids(file))
}

## ---- 2. Pairwise comparison ------------------------------------------

#' Percentage overlap between two closed intervals [s1,e1] and [s2,e2],
#' expressed as intersection / union * 100.
interval_overlap_pct <- function(s1, e1, s2, e2) {
  inter <- max(0, min(e1, e2) - max(s1, s2) + 1)
  union <- max(e1, e2) - min(s1, s2) + 1
  if (union <= 0) return(NA_real_)
  100 * inter / union
}

#' Compare two files' hit tables.
#'
#' @return list(n_shared = integer, pct_overlap = numeric or NA)
compare_hit_tables <- function(hits_a, hits_b) {
  shared_ids <- intersect(unique(hits_a$base_id), unique(hits_b$base_id))
  n_shared <- length(shared_ids)

  if (n_shared == 0) {
    return(list(n_shared = 0L, pct_overlap = NA_real_))
  }

  overlaps <- vapply(shared_ids, function(bid) {
    rows_a <- hits_a[hits_a$base_id == bid, , drop = FALSE]
    rows_b <- hits_b[hits_b$base_id == bid, , drop = FALSE]

    # If either side lacks coordinates for this ID, no overlap can be computed
    if (all(is.na(rows_a$start)) || all(is.na(rows_b$start))) {
      return(NA_real_)
    }

    # Multiple hits can share a base_id (e.g. multi-domain matches).
    # Take the best (maximum) overlap between any pair of intervals -
    # a natural "best match" pairing for that accession.
    best <- -Inf
    for (i in seq_len(nrow(rows_a))) {
      for (j in seq_len(nrow(rows_b))) {
        if (is.na(rows_a$start[i]) || is.na(rows_b$start[j])) next
        pct <- interval_overlap_pct(rows_a$start[i], rows_a$end[i],
                                     rows_b$start[j], rows_b$end[j])
        if (!is.na(pct) && pct > best) best <- pct
      }
    }
    if (is.infinite(best)) NA_real_ else best
  }, numeric(1))

  list(n_shared = n_shared, pct_overlap = mean(overlaps, na.rm = TRUE))
}

## ---- 3. Build matrices across a set of files --------------------------

#' Compute the two pairwise similarity matrices for a set of Stockholm files.
#'
#' @param files character vector of paths to Stockholm files
#' @return list(n_shared_matrix = matrix, pct_overlap_matrix = matrix)
compute_similarity_matrices <- function(files) {
  labels <- tools::file_path_sans_ext(basename(files))
  n <- length(files)

  hit_tables <- lapply(files, load_stockholm_hits)
  names(hit_tables) <- labels

  n_shared_matrix   <- matrix(NA_integer_, n, n, dimnames = list(labels, labels))
  pct_overlap_matrix <- matrix(NA_real_,    n, n, dimnames = list(labels, labels))

  for (i in seq_len(n)) {
    n_shared_matrix[i, i] <- nrow(hit_tables[[i]])
    pct_overlap_matrix[i, i] <- 100

    if (i < n) {
      for (j in (i + 1):n) {
        cmp <- compare_hit_tables(hit_tables[[i]], hit_tables[[j]])
        n_shared_matrix[i, j] <- n_shared_matrix[j, i] <- cmp$n_shared
        pct_overlap_matrix[i, j] <- pct_overlap_matrix[j, i] <- cmp$pct_overlap
      }
    }
  }

  list(n_shared_matrix = n_shared_matrix, pct_overlap_matrix = pct_overlap_matrix)
}

## ---- 4. Command-line entry point --------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript stockholm_pairwise_similarity.R <directory> [output_prefix]")
  }

  input_dir <- args[1]
  prefix    <- if (length(args) >= 2) args[2] else "stockholm_similarity"

  files <- list.files(input_dir, pattern = "\\.(sto|stk|stockholm)$",
                       full.names = TRUE, ignore.case = TRUE)
  if (length(files) < 2) {
    stop("Need at least two Stockholm files in: ", input_dir)
  }
  message("Found ", length(files), " Stockholm files.")

  result <- compute_similarity_matrices(files)

  write.csv(result$n_shared_matrix, paste0(prefix, "_n_shared.csv"))
  write.csv(result$pct_overlap_matrix, paste0(prefix, "_pct_overlap.csv"))

  message("Wrote ", prefix, "_n_shared.csv and ", prefix, "_pct_overlap.csv")
}

if (identical(environment(), globalenv()) && sys.nframe() == 0) {
  main()
}