% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rscrublet.R
\name{call_doublets}
\alias{call_doublets}
\title{Automatically sets a threshold for calling doublets,
but it's best to check this by running plot_doublet_histogram() afterwards and adjusting threshold by call_doublets with user defined \code{threshold} parameter}
\usage{
call_doublets(nnc, threshold = NULL, verbose = TRUE)
}
\arguments{
\item{nnc}{list, output of \link{scrub_doublets} function}

\item{threshold}{Doublet score threshold for calling a transcriptome
a doublet. If NULL, this is set automatically by looking
for the minimum between the two modes of the \code{nnc$doublet_scores_sim}
histogram. It is best practice to check the threshold visually
using the \code{nnc$doublet_scores_sim} histogram and/or based on
co-localization of predicted doublets in a 2-D embedding.}

\item{verbose}{If TRUE, print summary statistics.}
}
\value{
\code{nnc} list with additional items:
#' \itemize{
\item predicted_doublets - logical index of predicted doublet cells
\item z_scores_ - Z-score conveying confidence in doublet calls.
\item threshold - Doublet score threshold for calling a transcriptome a doublet.
\item detected_doublet_rate - Fraction of observed transcriptomes that have been called doublets.
\item detectable_doublet_fraction - Estimated fraction of doublets that are detectable, i.e., fraction of simulated doublets with doublet scores above \code{threshold}
\item overall_doublet_rate - Estimated overall doublet rate, \code{detected_doublet_rate/detectable_doublet_fraction}. Should agree (roughly) with \code{expected_doublet_rate}.
}
}
\description{
Automatically sets a threshold for calling doublets,
but it's best to check this by running plot_doublet_histogram() afterwards and adjusting threshold by call_doublets with user defined \code{threshold} parameter
}
\examples{
# run rscrublet of 8k pbmc example dataset
scrr = scrub_doublets(E_obs = pbmc8k,expected_doublet_rate=0.06,min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
# set threshould automatically 
scrr=call_doublets(scrr)
# examine score distribution
plot_doublet_histogram(scrr)
# find predicted doublets
rownames(pbmc8k)[scrr$predicted_doublets]
}
