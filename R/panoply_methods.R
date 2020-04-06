#
# panoply_methods.R
#
# declares exposed functions.
#
# 2019-12-17 WTR
#
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomics LLC 2019.                 ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics' express written consent.            ##
################################################################################
################################################################################

#' @import flowCore
#' @import flowFP
#' @import cluster
#' @import fields
#' @import Rtsne
#' @import KernSmooth
#' @import diptest
#' @importFrom grDevices chull
#' @importFrom graphics axis contour layout lines par pie points rect segments text
#' @importFrom methods as is new
#' @importFrom stats as.hclust cutree dist fivenum median quantile rnorm
#' @importFrom flowCore colnames identifier
NULL

#' @title panoply
#' @description This function wraps most of what's needed.  It computes a fingerprint
#' model, calculates the multivariate bin centers, imbeds them in a t-SNE map, and
#' performs clustering in t-SNE space.
#' @param fcs The data (either a flowFrame or a flowSet)
#' @param parameters The parameters in fcs used for analysis
#' @param nRecursions The number of recursions in calculating the fingerprint (default = 12)
#' @param perplexity The perplexity value used for the imbedding of the bins (default = 40)
#' @param nclust The number of clusters you want panoply to make.
#' @param merge Logical.  Merge categorically identical clusters.
#' @description Panoply implements a workflow of doing Cytometric Fingerprint (CF) binning of
#' a flowFrame (or flowSet)
#' using flowFP to a relatively high resolution (default nRecursions of 12 results
#' in 4096 bins).  The bin centroids are then computed, using the median value of all
#' of the events in the bin for each of the included parameters.  Next, these
#' multivariate bin centroids are imbedded in a 2-dimensional t-SNE map. In this case,
#' the dots in the t-SNE map correspond to CF bins, NOT individual events.  Finally,
#' The bins are clustered using aglommerative hierarchical clustering in the 2
#' t-SNE dimensions, with a user-specified number of clusters.
#'     The result is an object of class "panoply" that contains slots for all of
#' the important data elements.
#' @return An object of class panoply, with the following elements:
#' \describe{
#'   \item{mod}{The flowFPModel generated}
#'   \item{mfi}{A list representation of the bin centers}
#'   \item{centers}{The t-SNE map.  Dots represent bins}
#'   \item{map}{The t-SNE map.  Dots represent bins}
#'   \item{clustering}{A named list containing cluster membership of the bins}
#' }
#' @examples
#' # load the example data
#' load(system.file("extdata", "sampled_flowset_young.rda", package = "panoplyCF"))
#' pan_params = c(7:9, 11:22)
#' panoply(fs_young, parameters = pan_params, nclust = 15)
#' @export
panoply = function(fcs, parameters = NULL, nRecursions = 12, perplexity = 40, nclust = NULL, merge = TRUE) {
  # should not have to load libraries here, but did it to get vignette to knit
  # require(flowCore)
  # require(flowFP)
  if (is(fcs, "flowFrame")) {
    ff = fcs
  } else if (is(fcs, "flowSet")) {
    ff = suppressWarnings(as(fcs, "flowFrame"))
    flowCore::exprs(ff) = flowCore::exprs(ff)[,which(flowCore::colnames(flowCore::exprs(ff)) != "Original")]
  } else {
    stop("Argument fcs must either be a flowFrame or a flowSet\n")
  }
  # check parameters
  if (is.null(parameters)) {
    stop("Parameters must be either a numeric or character vector\n")
    if (is.numeric(parameters)) {
      parameters = flowCore::colnames(ff)[parameters]
    }
  }

  # check nclust
  if (is.null(nclust)) {
    stop("For now, you must specify how many clusters to generate\n")
  }

  message("computing fingerprint bins...")
  mod = flowFP::flowFPModel(ff, parameters = parameters, nRecursions = nRecursions)
  fp = flowFP::flowFP(ff, mod)
  message("calculating bin centers...")
  res = calculate_bin_phenotypes(fp = fp, fs = ff, method = "median")
  mfi = as.list(data.frame(t(res$center)))
  mat = t(res$center)

  # original t-SNE from PanoplyCF
  message("doing the t-SNE embedding...")

  map = do_tsne_reduction(mat, perplexity = perplexity, show = FALSE)   # higher perplexity seems to be a little nicer

  # cluster on the map
  message("agglomerative clustering in t-SNE space...")
  clst = cluster_map(map, k = nclust)

  # determining modality
  modality = parameter_modality(ff, parameters = parameters)

  panoply = list(mod = mod, mfi = mfi, centers = mat, map = map, clustering = clst, modality = modality)
  message("done.")

  # merging clusters
  if (merge) {
    panoply = merge_categorical_clusters(panoply = panoply, parameters = parameters)
    nclust = max(panoply$clustering$clst)
    message("merging clusters, resulting in ", nclust, " clusters...")
  }

  # calculate cluster centers, mostly for visualization
  c_centers = matrix (NA, nrow = 0, ncol = length(parameters))
  colnames(c_centers) = parameters
  for (i in 1:nclust) {
    idx = which(panoply$clustering$clst == i)
    c_centers = rbind(c_centers, distributions_bins(panoply, parameters, bin_indices = idx)$med)
  }
  panoply$clustering$c_centers = c_centers

  panoply = remap_clusters(panoply)

  class(panoply) = "panoply"
  return(panoply)
}

#' @title decorate_sample_panoply
#' @description This function makes a nice picture of the result of panoply.
#' @param fcs The data (either a flowFrame or a flowSet)
#' @param panoply_obj An object of type "panoply", the result of running panoply()
#' @param superclus TBD
#' @param colorscale logical.  Whether or not to draw the color wedge (default = TRUE)
#' @param superscale TBD
#' @param isLabeled TBD
#' @param cex Size of dots in peripheral panels.  Default 0.5
#' @param ... Additional graphical parameters (see par)
#' @description A picture of the result of running panoply().
#' @examples
#' decorate_sample_panoply(fs_young, pan)
#' @export
decorate_sample_panoply = function(fcs, panoply_obj, superclus=NULL,
                                   colorscale=TRUE, superscale=FALSE, isLabeled=FALSE,
                                   cex=0.5, ...) {
  if (is(fcs, "flowFrame")) {
    ff = fcs
  } else if (is(fcs, "flowSet")) {
    ff = suppressWarnings(as(fcs, "flowFrame"))
    flowCore::exprs(ff) = flowCore::exprs(ff)[,which(flowCore::colnames(flowCore::exprs(ff)) != "Original")]
  } else {
    stop("Argument fcs must either be a flowFrame or a flowSet\n")
  }

  mod = panoply_obj$mod
  map = panoply_obj$map
  mfi = panoply_obj$mfi
  clst = panoply_obj$clustering

  laymat = make_laymat(k = length(mfi), double = FALSE, allow_wedge = FALSE)
  layout(laymat)
  par(mar = c(0, 0, 0, 0) + 0.1)

  if (is.null(panoply_obj$clustering$func_phenotype)) {
    # make a cluster map
    count = NULL
    draw_cluster_map(map = map, clst = clst, superclus = superclus)
  } else {
    # show the functional cluster bubbles
    plot_tsne(panoply_obj, marker = "categorical",
              box = TRUE, cex = 50.0, proportional = TRUE, emph = TRUE,
              highlight_clusters = NULL, contours = TRUE, show_cluster_number = TRUE)
  }

  kde = bkde2D(map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
  min_value = 0
  max_value = 5
  par(mar = c(0, 0, 2, 0) + 0.1)
  mfi_lut = 1:length(mfi)
  for (i in 1:length(mfi)) {
    pname = names(mfi)[mfi_lut[i]]
    contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE,
            xaxt = 'n', yaxt = 'n',
            col = 'darkgray', main = pname, cex.main = 2)

    mod_min_value = min_value
    cols = pcolor(mfi[[mfi_lut[i]]], min_value = mod_min_value, max_value = max_value)
    points(map, cex = 0.9 * cex)
    points(map, col = cols, pch = 20, cex = cex)
  }
  if(colorscale){draw_color_scale(min_col_value = min_value, max_col_value = max_value)}
  if(superscale){
    if(length(superclus)<=12){
      RColorBrewer::display.brewer.pal(length(superclus),"Set3")
      if(isLabeled==TRUE){
        text(labels=names(superclus),x=1:length(superclus),y=1,cex=1.5,srt=90)
      }else{
        text(labels=1:length(superclus),x=1:length(superclus),y=1,cex=2.5)
      }
    }else{
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      cl<-col_vector[1:length(superclus)]
      pie(rep(1,length(superclus)), col=cl)
    }
  }
}

#' @title plot_panoply_spread
#' @description Draw a picture of the result of panoply for a given set of markers.
#' @param panoply The result of running panoply
#' @param markers Markers to include in the spread
#' @param mode Compute colors using either arithmetic or robust average
#' @param cex Scale factor for node size
#' @param proportional Logical.  Scale by the number of events in the cluster
#' @param emph Logical.  Emphasize each blob with a black line.
#' @param cex.lab Scale factor for titles of the individual markers
#' @param highlight_clusters Numeric vector indicating which clusters to highlight.
#' @param contours Logical. Plot the original t-SNE contours.
#' @param show_cluster_number Logical.  Print cluster numbers.
#' @return N/A.
#' @examples
#' plot_panoply_spread(panoply_obj)
#' @export
plot_panoply_spread = function(panoply, markers = NULL, mode = c("arithmetic", "robust"),
                             cex = 50.0, proportional = TRUE, emph = TRUE, cex.lab = 2,
                             highlight_clusters = NULL, contours = FALSE, show_cluster_number = FALSE) {
  plot_tsne_spread(panoply, markers, mode, cex, proportional, emph, highlight_clusters, contours, show_cluster_number)
  if (markers[1] != 'categorical') {
    draw_color_scale(cex.lab = cex.lab)
  } else {
    draw_cluster_legend(panoply_obj = panoply, cex.text = cex.lab)
  }
}


#' title map_functional_names
#' @description Based on a user-specified table, assign symbolic names to clusters
#' @param panoply_obj The result of running panoply()
#' @param defs_file A file containing the functional definitions.
#' @return A decorated panoply object
#' @export
map_functional_names = function(panoply_obj, defs_file) {
  fd = retrieve_categorical_definitions(defs_file)
  panoply_obj = assign_functional_names(panoply_obj, fd)
}

#' @title Map a Sample to a Panoply Model
#' @description This function determines, for a single sample, the number of cells in each cluster.
#' @param ff A sample flowFrame
#' @param panoply_obj An object of type "panoply", the result of running panoply()
#' @return  Per-cluster counts and fractions
#' @export
panoply_map_sample = function(ff, panoply_obj) {
  # apply the flowFPModel to the sample
  fp = flowFP(fcs = ff, model = panoply_obj$mod)

  # get the vector of event bin membership
  btag = tags(fp)[[1]]

  # assign event cluster membership
  nclust = max(panoply_obj$clustering$clst)
  nevents = nrow(ff)
  c_count = vector('numeric', length = nclust)
  for (i in 1:nclust) {
    bidx = panoply_obj$clustering$c_index[[i]]
    eidx = which(btag %in% bidx)
    c_count[i] = length(eidx)
  }

  # convert to percentages of total cells
  c_pctg = c_count / nevents

  # return the result
  return(list(counts = c_count, fractions = c_pctg))
}

#' @title Visualize Cluster Phenotypes
#' @description Draw a "phenobar" representation of a cluster phenotype.  Bars have
#' a height equal to the medial value of the parameter and are
#' color-coded.  Error flags represent first and third quartiles of the bin centers
#' belonging to the cluster.
#' @param panoply_obj An object of type "panoply", the result of running panoply()
#' @param parameters Which parameters to include in the plot (default = all parameters)
#' @param cluster Which cluster to plot.
#' @export
panoply_phenobars = function(panoply_obj,
                             parameters = colnames(panoply_obj$centers),
                             cluster = 1,
                             main = paste("Cluster", cluster)) {

  # make an empty plot
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '',
       main = main)
  wadeTools::ax(1, type = 'biexp')
  axis(side = 2, labels = parameters, at = 1:length(parameters), las = 1)

  centers = panoply_obj$centers

  # get the bin indices of the cluster
  idx = which(panoply_obj$clustering$clst == cluster)

  med_vec = vector(mode = 'numeric')
  q1_vec = vector(mode = 'numeric')
  q3_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    tmp = fivenum(centers[idx, i])
    med_vec[i] = tmp[3]
    q1_vec[i] = tmp[2]
    q3_vec[i] = tmp[4]
  }
  # draw the median
  col = pcolor(med_vec, min_value = 0, max_value = 5)
  add_bars(vals = med_vec, yvals = 1:length(parameters), col = col)


  # draw the flags
  for (i in 1:length(parameters)) {
    draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = NA, cex = 2, lwd = 2)

  }

}
