#! /usr/bin/env Rscript

################################################
# Score Inference Algorithm Submissions        #
# The Connectivity Map at The Broad Institute  #
# with HBS and Topcoder                        #
# Spring 2016                                  #
################################################

if (!(require('argparser'))) {
  message('argparser library not found, installing...')
  chooseCRANmirror(81)
  install.packages('argparser')
}

library(argparser)

# create arg parser for parsing command line arguments
parser <- arg_parser("
 Score inference algorithm results. Computes the Relative
 Inference Score (RIS) for a matrix of inferred gene expression
 values, and the average RIS across all genes. Outputs a summary table
 with the gene-level RIS values and a file, RIS.txt, with the
 mean RIS across all genes.

 example:
 cmap_scoring_funtion.R --inf_ds /path/to/inference/result.csv \
    --truth_ds /path/to/truth/matrix.csv \
    --reference_scores /path/to/reference/scores.csv \
    --out /desired/output/path
                     ")
parser <- add_argument(parser, "--inf_ds",
                       help="Path to matrix with inferred values.",
                       default=NA)
parser <- add_argument(parser, "--truth_ds",
                      help="Path to matrix with ground truth values.",
                       default=NA)
parser <- add_argument(parser, "--reference_scores",
                       help="Path to 1-column CSV containing a vector of gene-level
                         scores derived on the existing CMap inference model. If this is
                         not provided, no relative scoring is done.",
                       default=NA)
parser <- add_argument(parser, "--gene_dim",
                       help="Indicates which dimension in --inf_ds and --truth_ds 
                         corresponds to genes. Assumed to be the same for both matrices.
                         Must be one of ['row', 'column'].",
                       default="row")
parser <- add_argument(parser, "--headerless",
                       help="Indicates whether the matrices are without row and column
                         identifiers. If not, it is assumed that the rows and columns
                         of --inf_ds and --truth_ds are in sync.",
                       flag=T)
parser <- add_argument(parser, "--make_plots",
                       help="Whether to generate plots.",
                       flag=T)
parser <- add_argument(parser, "--out", help="output path", default=getwd())
args <- parse_args(parser, argv=commandArgs(trailingOnly=T))

# if args are valid, load more libraries
msg.trap <- capture.output(suppressMessages(library(roller)))

### Function Definitions ###

# fast correlations
fast_cor <- function(X, Y) {
  # convert to ranks by column
  X <- apply(X, 2, rank, ties.method = "average")
  Y <- apply(Y, 2, rank, ties.method = "average")
  # zscore columns
  X <- scale(X)
  Y <- scale(Y)
  # spearman
  nsample <- nrow(X)
  corr <- (t(X) %*% Y) / (nsample - 1)
  return(corr)
}

# THE RELATIVE SCORING FUNCTION
# computes the relative imputation score (RIS)
get_relative_score <- function(ref_scores, test_scores, coefficient=1e6) {
  numerator <- 2 - ref_scores
  denominator <- 2 - test_scores
  scores <- numerator / denominator
  scores <- coefficient * scores
  return(scores)
}

# rescale the relative score
# assumes current_score is the aggregated RIS
get_rescaled_score <- function(current_score, score_ceiling, par_score=1e6) {
  # ScoreScaled = CurrentScore (ScoringCeiling - 1M)/(ScoringCeiling - CurrentScore)
  score_scaled <- current_score * ((score_ceiling - par_score) / (score_ceiling - current_score))
  return(score_scaled)
}

### Main Program ###

# make sure required args are supplied
req_args <- c("inf_ds", "truth_ds", "gene_dim")
missing_args <- setdiff(req_args, names(args))
if (length(missing_args) > 0) {
  stop(paste("the following required arguments are missing:\n",
             paste(missing_args, collapse="\n")))
}

# make sure the gene_dim argument is valid
if (!args$gene_dim %in% c("row", "column")) {
  stop("gene_dim must be either 'row' or 'column'")
}

# set up output path
outpath <- args$out
if (!file.exists(outpath)) dir.create(outpath)

# read the user's submission
# expected to be a CSV containing
# the predicted expression levels.
# also read the ground truth data
# expected to be a CSV containing
# the ground truth expression levels
message("reading data...")
if (args$headerless) {
  inf <- as.matrix(read.csv(args$inf_ds, header=F))
  truth <- as.matrix(read.csv(args$truth_ds, header=F))
} else {
  inf <- as.matrix(read.csv(args$inf_ds, row.names=1))
  truth <- as.matrix(read.csv(args$truth_ds, row.names=1))
  # make sure rows and samples are common
  different_rows <- setdiff(rownames(inf), rownames(truth))
  different_cols <- setdiff(colnames(inf), colnames(truth)) 
  if (length(different_rows) > 0) {
    stop("matrices have different row ids")
  }
  if (length(different_cols) > 0) {
    stop("matrices have different column ids")
  }
  inf <- inf[rownames(truth), colnames(truth)]
}
message("done")

# make sure matrices are the same dimensions
if (nrow(inf) != nrow(truth) | ncol(inf) != ncol(truth)) {
  stop("matrices are of different dimensions")
}

# need to get genes on the columns, so tranpose
# matrices if necessary
if (args$gene_dim == "row") {
  message("transposing matrices so genes are along columns")
  inf <- t(inf)
  truth <- t(truth)
}

# compute correlations
message("computing correlations...")
corr <- fast_cor(inf, truth)
message("done")

# rank across columns (inf genes on rows)
# asks how well the inferred gene pulls up its
# measured counterpart
# higher values get higher rank
rank_by_row <- apply(corr, 1, function(x) rank(x))

# convert to fraction ranks
frac_ranks <- rank_by_row / nrow(corr)

# collate the absolute correlation and the
# fraction rank into a single table
df <- data.frame(correlation = diag(corr),
                 fraction_rank = diag(frac_ranks))
if (!args$headerless & !is.null(names(diag(corr)))) {
  df$gene <- names(diag(corr))
}

# average the correlation and fraction rank to get the score
# threshold negative correlations to zero
df$thresholded_correlation <- sapply(df$correlation, function(x) max(x, 0))
df$score_test <- (df$thresholded_correlation + df$fraction_rank) / 2 
mean_score <- mean(df$score_test, na.rm=T)

# start printing a summary to the screen
message("---------- Summary ----------")
message(paste("Mean Score:", round(mean_score, 2)))

# reference_scores provided
# get the reference score relative to baseline
if (!is.na(args$reference_scores)) {
  # read the baseline performance table
  # if --headerless, expect this to be a single-column table
  if (args$headerless) {
    base_performance <- as.vector(read.csv(args$reference_scores, header=F)[, 1])
    df$score_ref <- base_performance
  } else {
    base_performance <- read.csv(args$reference_scores)
    df <- merge(df, base_performance, by="gene")
    setnames(df, "score", "score_ref")
  }
  
  df$RIS <- get_relative_score(df$score_ref, df$score_test)

  # figure out the max possible score by supplying
  # a vector of all perfect scores (all 1)
  best_ris <- get_relative_score(df$score_ref, rep(1, nrow(df)))
  mean_best_ris <- mean(best_ris)

  # compute some summary stats
  mean_RIS <- mean(df$RIS, na.rm=T)

  # rescaled mean_RIS
  rescaled_mean_RIS <- get_rescaled_score(mean_RIS, mean_best_ris)

  # dump a file with the mean RIS
  write(paste("RIS:", round(rescaled_mean_RIS, 2)), file=paste(outpath, "RIS.txt", sep="/"),
        append = T)

  # keep printing summary stats
  message(paste("Mean RIS:", round(mean_RIS, 2)))
  message(paste("Theoretical Best RIS:", round(mean_best_ris, 2)))
  message(paste("Mean RIS Rescaled:", round(rescaled_mean_RIS, 2)))
}

# dump the summary table
write.csv(df, paste(outpath, "summary_table.csv", sep="/"),
          row.names=F, quote=F)

# were we asked to make plots?
if (args$make_plots) {
  
  # figure out the plot title based on whether we computed
  # relative scores
  if (exists("mean_RIS")) {
    title <- paste(
      "Mean Correlation:", round(mean(df$thresholded_correlation, na.rm=T), 2), "\n",
      "Mean Fraction Rank:", round(mean(df$fraction_rank, na.rm=T), 2), "\n",
      "Mean RIS:", round(mean_RIS, 2), "\n",
      "Mean RIS Rescaled:", round(rescaled_mean_RIS, 2)
    )
  } else {
    title <- paste(
      "Mean Correlation:", round(mean(df$thresholded_correlation, na.rm=T), 2), "\n",
      "Mean Fraction Rank:", round(mean(df$fraction_rank, na.rm=T), 2), "\n"
    )
  }
  
  # correlation vs. fraction rank
  # scatter w/ marginal histograms
  p <- ggplot(df, aes(x=thresholded_correlation, y=fraction_rank))
  p <- p + theme_classic(base_size = 30)
  p <- p + geom_point(color="black", shape=1)
  p <- p + xlim(0, 1) + ylim(0, 1)
  p <- p + xlab("Correlation") + ylab("Fraction Rank")
  p <- p + theme(panel.border=element_rect(fill=NA))
  p <- p + ggtitle(title)
  png(path.join(outpath, "rank_vs_correlation_scatter.png"), height=768, width=768, units="px")
  # marginal histograms
  print(ggExtra::ggMarginal(p,
                            type = "histogram",
                            color = "black",
                            fill = "gray"))
  dev.off()
  
  # their gene-level scores vs. reference
  # scatter w/ marginal histograms
  # (if we have reference scores)
  if (!is.null(df$score_ref)) {
    p <- ggplot(df, aes(x=score_ref, y=score_test))
    p <- p + theme_classic(base_size = 30)
    p <- p + geom_point(color="black", shape=1)
    p <- p + xlim(0, 1) + ylim(0, 1)
    p <- p + xlab("Reference Score") + ylab("Test Score")
    p <- p + geom_abline(slope=1, intercept=0, linetype="dashed", color="red")
    p <- p + theme(panel.border=element_rect(fill=NA))
    p <- p + ggtitle(title)
    png(path.join(outpath, "test_vs_ref_scatter.png"), height=768, width=768, units="px")
    # marginal histograms
    print(ggExtra::ggMarginal(p,
                              type = "histogram",
                              color = "black",
                              fill = "gray"))
    dev.off()
  }

}
