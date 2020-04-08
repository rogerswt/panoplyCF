#
# panoply_utils.R
#
# These are functions to support PanoplyCF.  Not intended to be exposed to the
# user.
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

# input is a flowFP (fp) and the corresponding flowSet (fs)
#  method = "median" returns median +- quartiles.
#  method = "mean" returns mean +- standard deviation
calculate_bin_phenotypes = function(fp, fs, method=c("median", "mean")) {
  parameters = parameters(fp)
  n_bins = 2 ^ nRecursions(fp)
  n_parameters = length(parameters)
  center = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(center) = parameters
  range = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(range) = parameters

  # for convenience, lump the frames in fs together and recalculate the fingerprint
  ff = as(fs, "flowFrame")
  fp = flowFP(fcs = ff, model = as(fp, "flowFPModel"))
  for (i in 1:n_bins) {
    idx = which(tags(fp)[[1]] == i)
    if (length(idx) == 0) {
      next
    }

    for (j in 1:n_parameters) {
      p = parameters[j]
      vals = exprs(ff)[idx, p]
      if (method == "median") {
        center[j, i] = median(vals, na.rm = TRUE)
        range[j, i] = (quantile(x = vals, probs = 0.75, na.rm = TRUE) -
                         quantile(x = vals, probs = 0.25, na.rm = TRUE)) / 2
      } else {

      }
    }
  }
  return(list(center = center, range = range))
}

# given a collection of bin centers, perform T-SNE dimensionality reduction
do_tsne_reduction = function(centers, perplexity = 30, show=FALSE) {
  set.seed(137)   # so we'll get the same map for the same data
  res = Rtsne(dist(centers), perplexity = perplexity)$Y
  colnames(res) = c("t_sne_1", "t_sne_2")
  if (show) {
    plot(res, pch = 20, col = 'red')
    points(res)
  }

  res

}

# map a sample to our Panoply model
map_sample_panoply = function(tubes, mod, map
                              , min_val = -5, max_val = -1
                              , norm_fac = 1.0
                              , x_leg = -25, y_leg = -25
                              , ...) {
  kde = bkde2D(map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))

  # aggregate the 3 tubes, and extract the fp parameters
  if (is.list(tubes)) {
    flist = extract_common(tubes, parameters(mod))
    fs = as(flist, "flowSet")
    ff = as(fs, "flowFrame")
  } else {
    ff = tubes
  }
  n_tot = nrow(ff)
  tmp = flowFP(ff, mod)
  count = counts(tmp) * norm_fac
  cnt = as.vector(log10((count + 1) / n_tot))
  # draw the map using a saturation ramp (white is low, dark is high)
  # spread values from min_val -- max_val (white -- dark)
  # anything less than -5 is white, anything above
  n_color = 100
  lut_col = heat.color.sat()(n_color)
  idx_col = round(n_color * (cnt - min_val) / (max_val - min_val))
  idx_col[which(idx_col < 1)] = 1
  idx_col[which(idx_col > n_color)] = n_color
  cols = lut_col[idx_col]
  if (is.null(kde)) {
    plot(map, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', ...)
  } else {
    contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE,
            xaxt = 'n', yaxt = 'n',
            col = 'darkgray', ...)
    points(map, ...)

  }
  points(map, pch = 20, col = cols, ...)

  insert_wedge(x_ll = x_leg, y_ll = y_leg, width = 2, height = 15, colvec = lut_col,
               min_val = min_val, max_val = max_val)
  invisible(count)
}

make_laymat = function(k, double = FALSE, allow_wedge = FALSE) {
  #cat(paste("asking for k of ",k,"\n"))
  if (k == 4) {    # special case
    laymat = matrix(0, nrow = 3, ncol = 3)
    laymat[2, 2] = 1
    laymat[1, 1] = 2
    laymat[1, 3] = 3
    laymat[3, 1] = 4
    laymat[3, 3] = 5
    laymat[3, 2] = 6

    return(laymat)
  }
  if (double) {
    n = ceiling(k/8 + 2)
    if (n <= 5) {n = 6}
  } else {
    n = ceiling(k/4 + 1)
  }

  laymat = matrix(0, nrow = n, ncol = n)

  # make the central figure
  if (double) {
    for (i in 3:(n - 2)) {
      for (j in 3:(n - 2)) {
        laymat[i, j] = 1
      }
    }
  } else {
    for (i in 2:(n - 1)) {
      for (j in 2:(n - 1)) {
        laymat[i, j] = 1
      }
    }
  }

  # top
  if (double) {
    laymat[1, ] = (1:n)  + 1
    laymat[2, ] = ((n + 1):(2 * n)) + 1
  } else {
    laymat[1, ] = (1:n)  + 1
  }

  # middle
  if (double) {
    m = (2 * n + 1) + 1
    for (i in 3:(n - 2)) {
      for (j in c(1, 2, n - 1, n)) {
        laymat[i, j] = m
        m = m + 1
      }
    }
  } else {
    if(k==4){
      laymat[2, 2]<-4
      return(laymat)
    }else{
      m = (n + 1) + 1
      for (i in 2:(n - 1)) {
        for (j in c(1, n)) {
          laymat[i, j] = m
          m = m + 1
        }
      }
    }
  }

  # bottom
  if (double) {
    laymat[(n - 1), ] = m:(m + n - 1)
    laymat[n, ]       = (m + n):(m + (2 * n) - 1)
  } else {
    laymat[n, ] = m:(m + n - 1)
  }

  if (allow_wedge) {
    if (max(laymat) == k + 1) {
      cat("warning, not enough room for wedge\n")
    }
    laymat[laymat > k + 1] = 0
    laymat[n, n] = k + 2
  }

  laymat
}

# assumes biexp vert scale
draw_color_scale = function(min_col_value = 0, max_col_value = 5, ...) {
  ll = -0.5
  ul = bx(262143)

  vec = seq(ll, ul, length.out = 500)
  cols = pcolor(pvalue = vec, min_value = min_col_value, max_value = max_col_value)

  opar = par(mar = c(0, 5, 0, 0) + .1)
  plot(0, 0, pch = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "Fluorescence Intensity",
       xlim = c(0, 5), ylim = c(ll, ul), ...
  )
  for (i in 1:length(vec)) {
    y = vec[i]
    segments(x0 = 0, y0 = y, x1 = 1, y1 = y, col = cols[i], lwd = 3)
  }
  ax(axis = 2, instrument = 'diva', type = 'biexp', ...)
  par(opar)
}

# draw legend for categorical labeling
draw_cluster_legend = function(panoply_obj, cex.text = 1.25) {
  opar = par(mar = c(.1, .1, .1, .1))
  plot(0, 0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(0, 1), ylim = c(0, 1), bty = 'n')
  xcol = 0.05
  wcol = .1
  xtext = 0.18
  subsets = rownames(panoply_obj$labels$definitions)
  cols = panoply_obj$labels$colors

  n_items = length(subsets)
  y = seq(.05, .95, length.out = n_items)
  ht = .8 * .5 / n_items
  for (i in 1:n_items) {
    rect(xleft = xcol, ybottom = y[i] - ht, xright = xcol + wcol, ytop = y[i] + ht, col = cols[i])
    text(x = xtext, y = y[i], labels = subsets[i], pos = 4, cex = cex.text)
  }
  par(opar)
}

cluster_map = function(map, h = NULL, k = NULL) {
  if (is.null(h) & is.null(k)) {
    stop("Must provide EITHER n or k\n")
  }
  ag = agnes(map)
  clst = cutree(as.hclust(ag), h = h, k = k)
  n_clust = max(clst)

  centers = matrix(NA, nrow = n_clust, ncol = ncol(map))
  boundaries = list()
  cindex = list()
  colnames(centers) = c("t_sne_1", "t_sne_2")
  for (i in 1:n_clust) {
    idx = which(clst == i)
    cindex[[i]] = idx
    centers[i, ] = c(median(map[idx, 1]), median(map[idx, 2]))
    boundaries[[i]] = get.hull(map[idx, ])
  }

  return(list(clst = clst, c_index = cindex, centers = centers, boundaries = boundaries))
}

draw_cluster_map = function(map, clst, dot_col='gray', superclus=NULL) {
  kde = bkde2D(map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
  contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE, col = 'darkgray', add = FALSE, xaxt = 'n', yaxt = 'n')
  if (!is.null(superclus)) {
    # color-code the superclusters
    #cl = hsv(h = seq(0, .6667, length.out = length(superclus)), s = 1, v = .75)
    if(length(superclus)<=12){
      cl <- RColorBrewer::brewer.pal(length(superclus),"Set3")
    }else{
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      cl<-col_vector[1:length(superclus)]
    }
    dot_col = rep('gray', length.out = nrow(map))
    for (i in 1:length(superclus)) {
      idx = bins_in_supercluster(clst = clst, sc = superclus[[i]])
      dot_col[idx] = cl[i]
    }
  }
  points(map, pch = 20, col = dot_col, cex = 2)
  for (i in 1:nrow(clst$centers)) {
    lines(clst$boundaries[[i]], col = 'darkgreen')
  }
  text(x = clst$centers[, 1], y = clst$centers[, 2], labels = 1:nrow(clst$centers), vfont = c("serif", "bold"), cex = 1.5)

}

# insert a color scale wedge inside a plot
insert_wedge = function(x_ll, y_ll, width, height, colvec, min_val, max_val) {

  n_color = length(colvec)
  x0 = x_ll
  x1 = x_ll + width
  y = y_ll
  delta_y = height / n_color
  for (i in 1:n_color) {
    segments(x0, y, x1, y, col = colvec[i], lwd = 5)
    y = y + delta_y
  }
  y_lab = pretty(seq(min_val, max_val, length.out = n_color))
  x_lab = x_ll - .1 * width
  intrvl = height / (length(y_lab) - 1)
  y = y_ll
  for (i in 1:length(y_lab)) {
    text(x = x_lab, y = y, labels = y_lab[i], pos = 2, cex = 0.75)
    y = y + intrvl
  }
  rect(xleft = x_ll, ybottom = y_ll, xright = x_ll + width, ytop = y_ll + height, border = 'black')

}

# define a color function to transform a parameter value
pcolor = function(pvalue, min_value = 0, max_value = 4) {
  top_col = 'red'
  bot_col = 'darkgreen'
  mid_col = 'yellow'
  zero_col = 'darkgray'
  zero_col = bot_col

  len = length(pvalue)
  if (length(which(pvalue < min_value)) == len) {
    col_values = rep(zero_col, len)
  } else if (length(which(pvalue > max_value)) == len) {
    col_values = rep(top_col, len)
  } else {
    pvalue[pvalue < min_value] = min_value
    pvalue[pvalue > max_value] = max_value

    col_values = color.scale(
      z = pvalue, col = two.colors(
        n = 100,
        start = bot_col,
        end = top_col,
        middle = mid_col),
      zlim = c(min_value, max_value)
    )
    col_values[which(pvalue <= min_value)] = zero_col
  }


  col_values
}

# do some parallel coordinate plots
parallel_pheno = function(mfi, idx_bin = 1:length(mfi[[1]]), parameters = names(mfi), col = 'red', mfi_colors = FALSE, show_contours = TRUE, bars = TRUE, ...) {

  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  if (show_contours) {
    cont_col = 'darkgray'
  } else {
    cont_col = 'white'
  }
  mfi = data.frame(mfi)
  # hyphens in parameter names get turned into dots here.  Change back
  colnames(mfi) = sub(pattern = ".", replacement = "-", x = colnames(mfi), fixed = TRUE)

  parallel_contours(mfi, parameters = parameters, col = cont_col, ...)
  mfi = mfi[idx_bin, parameters]

  med_vec = vector(mode = 'numeric')
  q1_vec = vector(mode = 'numeric')
  q3_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    tmp = fivenum(mfi[, i])
    med_vec[i] = tmp[3]
    q1_vec[i] = tmp[2]
    q3_vec[i] = tmp[4]
  }
  # draw the median
  if (bars) {
    if (mfi_colors) {
      col = pcolor(med_vec, min_value = 0, max_value = 5)
    }
    add_bars(vals = med_vec, yvals = 1:length(parameters), col = col)
  } else {
    x = med_vec
    y = 1:length(parameters)
    lines(x, y, col = col, lwd = 3)
  }

  # draw the flags
  for (i in 1:length(parameters)) {
    if (bars) {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = NA, cex = 2, lwd = 2)
    } else {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = med_vec[i], cex = 2, lwd = 2)
    }
  }
}

parallel_contours = function(mfi, parameters = names(mfi), col = 'blue', ...) {
  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '')
  ax(1, instrument = 'diva', type = 'biexp')
  axis(side = 2, labels = pnames, at = 1:length(pnames), las = 1, ...)

  x = matrix(NA, nrow = nrow(mfi) * length(parameters), ncol = 2)
  k = 1
  for (i in 1:nrow(mfi)) {
    for (p in 1:length(parameters)) {
      x[k, 1] = mfi[i, p]
      x[k, 2] = p
      k = k + 1
    }
  }
  kde = bkde2D(x = x, bandwidth = c(.1, 1), gridsize = c(501, 501))
  kde$fhat = kde$fhat / max(kde$fhat)
  contour(x = kde$x1, y = kde$x2, z = kde$fhat, col = col,
          xaxt = 'n', yaxt = 'n',
          drawlabels = FALSE,
          # levels = seq(.01, .2, length.out = 20),
          add = TRUE)

}

draw_y_grid = function(lo = -4, hi = 0) {
  minor = seq(2, 9, by = 1)
  major = 10 ^ seq(lo, hi, by = 1)
  for (maj in major) {
    yline(maj, col = "darkgray", lwd = 2)
    for (mn in minor) {
      yline(maj * mn, col = 'lightgray', lwd = 1)
    }
  }
}

add_bars = function(vals, yvals, col) {
  hw = 0.4
  if (length(col) == 1) {
    col = rep(col, length(vals))
  }
  for (i in 1:length(vals)) {
    rect(xleft = 0, ybottom = yvals[i] - hw, xright = vals[i], ytop = yvals[i] + hw, col = col[i], border = 'black')
  }
}

draw_flag = function(y, q1, q3, med = NA, ...) {
  segments(x0 = q1, y0 = y, x1 = q3, y1 = y, ...)
  if (!is.na(med)) {points(med, y, pch = 20, ...)}
}

get.hull <- function (blob) {
  # 2018-08-27 WTR - handle edge case that blob is a single point
  if (is.vector(blob)) {
    cnames= names(blob)
    blob = matrix(blob, nrow = 1, ncol = 2)
    colnames(blob) = cnames
    return(blob)
  }
  if(!is.matrix(blob)) {
    cnames = names(blob)
    blob = matrix(blob, ncol = 2)
    colnames(blob) = cnames
    return(blob)
  }
  x <- blob[,1]
  y <- blob[,2]
  hull <- chull(x, y)
  poly <- matrix(c(x[hull], y[hull]), nrow=length(hull), ncol=2)
  # get rid of dupes
  poly <- unique(poly)
  # close the contour
  poly <- rbind(poly, poly[1,])
  colnames(poly) <- colnames(blob)
  return(poly)
}

# calculate statistics for a collection of bins
distributions_bins = function(panoply_obj, parameters = colnames(panoply_obj$centers), bin_indices) {
  centers = panoply_obj$centers
  idx = bin_indices

  med_vec = vector(mode = 'numeric')
  q1_vec = vector(mode = 'numeric')
  q3_vec = vector(mode = 'numeric')
  mn_vec = lo_vec = hi_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    tmp = fivenum(centers[idx,parameters[i]])
    med_vec[i] = tmp[3]
    q1_vec[i] = tmp[2]
    q3_vec[i] = tmp[4]
    mn_vec[i] = mean(centers[idx,parameters[i]])
    sdev = sd(centers[idx,parameters[i]])
    lo_vec[i] = mn_vec[i] - 0.5 * sdev
    hi_vec[i] = mn_vec[i] + 0.5 * sdev
  }
  names(med_vec) = names(q1_vec) = names(q3_vec) = parameters

  return(list(med = med_vec, q1 = q1_vec, q3 = q3_vec, mn = mn_vec, lo = lo_vec, hi = hi_vec))
}

# find the location of each cluster within the original map
remap_clusters = function(panoply_obj) {
  n_clust = length(panoply_obj$clustering$c_index)

  remap = matrix(NA, nrow = n_clust, ncol = 2)
  colnames(remap) = c("t_sne_1", "t_sne_2")

  for (i in 1:n_clust) {
    # get the bins in the cluster
    idx = panoply_obj$clustering$c_index[[i]]
    remap[i, ] = c(median(panoply_obj$map[idx, 1]), median(panoply_obj$map[idx, 2]))
  }

  panoply_obj$clustering$remap = remap

  panoply_obj
}

# label a cluster categorically
# If the median for the cluster is above the parameter threshold, label it hi, otherwise lo
categorical_phenotype = function(panoply_obj, parameters = colnames(panoply_obj$centers), cluster = 1) {
  idx = panoply_obj$clustering$c_index[[cluster]]
  res = distributions_bins(panoply_obj, bin_indices = idx)
  # label each parameter as either lo or hi
  p_category = rep("lo", length = length(parameters))
  names(p_category) = parameters
  for (i in 1:length(parameters)) {
    if (res$mn[i] > panoply_obj$modality$thresh[i]) {p_category[i] = "hi"}
  }
  p_category = factor(p_category, levels = c("lo", "hi"))

  p_category
}

compare_categories = function(c1, c2) {
  n_param = length(c1)
  if (length(which(c1 == c2)) == n_param) {
    res = TRUE
  } else {
    res = FALSE
  }
  res
}

merge_categorical_clusters = function(panoply_obj, parameters = colnames(panoply_obj$centers)) {
  # get the categorical mapping
  n_clust = max(panoply_obj$clustering$clst)
  categ = list()
  for (i in 1:n_clust) {
    categ[[i]] = categorical_phenotype(panoply_obj = panoply_obj, parameters = parameters, cluster = i)
  }

  # roll through and create clusters of clusters
  cmerge = list()
  phenotype = list()
  # cvec is a vector of cluster indices.  When a cluster joins a merge, it's removed from this vector
  cvec = 1:n_clust
  k = 1

  while (length(cvec) > 1) {
    # get the head of the list of remaining clusters
    ith = cvec[1]
    cmerge[[k]] = ith                          # add ith to the next merge
    phenotype[[k]] = categ[[ith]]              # record the phenotype
    cvec = cvec[which(cvec != ith)]            # remove it from cvec
    j = 1
    while (j <= length(cvec)) {
      jth = cvec[j]
      if (compare_categories(categ[[ith]], categ[[jth]])) {
        cmerge[[k]] = append(cmerge[[k]], jth)     # add jth cluster to cmerge
        cvec = cvec[which(cvec != jth)]            # remove jth cluster from cvec
        j = j - 1                                  # don't skip next element
      }
      j = j + 1
    }
    k = k + 1
  }
  # handle the last cluster
  if (length(cvec) > 0) {
    cmerge[[k]] = cvec
    phenotype[[k]] = categ[[cvec]]
  }

  # replace clustering slot with the merged result
  orig_clustering = panoply_obj$clustering
  c_index = list()
  for (i in 1:length(cmerge)) {
    c_index[[i]] = vector(mode = 'numeric')
    for (j in 1:length(cmerge[[i]])) {
      c_index[[i]] = append(c_index[[i]], orig_clustering$c_index[[cmerge[[i]][j]]])
    }
  }
  nbins = 2 ^ nRecursions(panoply_obj$mod)
  clst = rep(NA, length = nbins)
  for (i in 1:length(c_index)) {
    clst[c_index[[i]]] = i
  }

  n_clust = length(cmerge)
  map = panoply_obj$map
  centers = matrix(NA, nrow = n_clust, ncol = ncol(map))
  boundaries = list()
  cindex = list()
  colnames(centers) = c("t_sne_1", "t_sne_2")
  for (i in 1:n_clust) {
    idx = which(clst == i)
    cindex[[i]] = idx
    centers[i, ] = c(median(map[idx, 1]), median(map[idx, 2]))
    boundaries[[i]] = get.hull(map[idx, ])
  }

  clustering = list(clst = clst, c_index = c_index, centers = centers, boundaries = boundaries, phenotype = phenotype, func_phenotype = NULL)
  panoply_obj$clustering = clustering

  panoply_obj
}

# use the Hartigan dip-test to estimate modality for each parameter
# calculate thresholds to be used to determine lo/hi categorization
parameter_modality = function(ff, parameters = detect_fl_parameters(ff), crit = .2) {
  # down-sample ff due to dip.test limitation
  if (is(ff, "flowSet")) {ff = as(ff, "flowFrame")}
  # Next two lines were to avoid a max_n message in dip.test, but sampling
  # leads to non-deterministic results for cluster merging.
  # max_n = 71999
  # if(nrow(ff) > max_n) {ff = Subset(ff, sampleFilter(max_n))}

  pv = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    pv[i] = suppressMessages(dip.test(exprs(ff)[, parameters[i]])$p.value)
  }

  names(pv) = parameters
  unimodal = pv > crit

  # calculate hi/lo thresholds, conditioned on modality
  thresh = vector(mode = 'numeric')
  # multimodal thresholds
  for (i in which(!unimodal)) {
    kde = bkde(exprs(ff)[, parameters[i]])
    res = find.local.minima(kde)
    # find the lowest one above bx(1000)
    thresh[i] = min(res$x)
    if (length(res$x) > 1) {
      thresh[i] = min(res$x[which(res$x > bx(1000))])
    }
  }
  # unimodal thresholds defined as greater than mean + 1sd
  for (i in which(unimodal)) {
    x = exprs(ff)[, parameters[i]]
    mn = mean(x)
    sdev = sd(x)
    thresh[i] = mn + sdev
  }
  names(unimodal) = parameters
  names(thresh) = parameters

  return(list(unimodal = unimodal, thresh = thresh, p.value = pv))
}

# Use a spreadsheet to define functional phenotypes.
# Rows are functional phenotypes (for example, CD4_EM)
# columns are markers (e.g. CD4, CD8, CD3, ...)
# Entries are either 'lo', 'hi', or 'dc' (indicating don't care)
# rownames will define the subset names
# colnames MUST match the marker labeling in the data set
# optional column labeled "Colors" will contain color coding of subsets
retrieve_categorical_definitions = function(file) {
  tab = read.csv(file, row.names = 1, as.is = TRUE)

  # now turn the columns into properly ordered factors
  idx = which(tolower(colnames(tab)) != "color")
  icol = (1:ncol(tab))[-idx]
  # for (i in idx) {
  #   tab[, i] = factor(tab[, i], levels = c("lo", "hi", "dc"))
  if (length(icol) == 1) {
    colors = tab[, icol]
  } else {
    colors = rep(NA, nrow(tab))
  }
  # }

  return(list(defs = tab[, idx], colors = colors))
}

assign_functional_names = function(panoply_obj, functional_definitions) {
  fd = as.matrix(functional_definitions$defs)
  fc = functional_definitions$colors
  # add to the panoply object
  panoply_obj$labels$definitions = fd
  panoply_obj$labels$colors = fc

  n_clust = length(panoply_obj$clustering$c_index)
  func_phenotype = rep('unassigned', length = n_clust)
  func_color = rep("black", length = n_clust)
  pheno = panoply_obj$clustering$phenotype
  for (i in 1:n_clust) {
    for (j in 1:nrow(fd)) {
      idx = which(fd[j, ] != "dc")
      cmarkers = colnames(fd)[idx]
      def = fd[j, idx]
      # extract values for cluster
      ctype = as.character(pheno[[i]][cmarkers])
      if (length(which(ctype == def)) == length(idx)) {
        func_phenotype[i] = rownames(fd)[j]
        func_color[i] = fc[j]
        break
      }
    }
  }

  # tack onto panoply_obj
  panoply_obj$clustering$func_phenotype = func_phenotype
  panoply_obj$clustering$func_color = func_color

  panoply_obj
}

plot_tsne = function(panoply_obj, marker, mode = c("arithmetic", "robust"),
                     box = TRUE, cex = 50.0, proportional = TRUE, emph = TRUE,
                     highlight_clusters = NULL, contours = FALSE, show_cluster_number = NULL) {
  mode = match.arg(mode)
  centers = panoply_obj$centers
  n_clust = length(panoply_obj$clustering$c_index)
  if (marker == "categorical") {
    cols = panoply_obj$clustering$func_color
  } else {
    cols = pcolor(panoply_obj$clustering$c_centers[, marker], min_value = 0.0, max_value = 5.0)
  }

  map = panoply_obj$clustering$remap
  if (proportional) {
    size = vector('numeric')
    for (i in 1:n_clust) {
      size[i] = length(panoply_obj$clustering$c_index[[i]])
    }
    size = sqrt(size / (2 ^ nRecursions(panoply_obj$mod)))
    cex = cex * size
    cexbg = 1.05 * cex
  }

  if (!is.null(highlight_clusters)) {
    hcol = rep('black', length = n_clust)
    hcol[highlight_clusters] = 'magenta'
    cexbg[highlight_clusters] = 2.0 * cexbg[highlight_clusters]
  } else {
    hcol = rep('black', length = n_clust)
  }

  # plot largest first
  srt = sort(size, decreasing = TRUE, index.return = TRUE)$ix
  map = map[srt, ]
  cols = cols[srt]
  cex = cex[srt]
  cexbg = cexbg[srt]
  hcol = hcol[srt]

  bty = ifelse(box, 'o', 'n')
  tit = ifelse(marker == 'categorical', '', marker)
  if (emph) {
    xlim = c(min(map[, 1]), max(map[, 1]))
    ylim = c(min(map[, 2]), max(map[, 2]))
    plot(0, 0, pch = '', , bty = bty, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = xlim, ylim = ylim, main = tit)
    if (contours) {
      kde = bkde2D(panoply_obj$map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
      contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE, col = 'darkgray', add = TRUE, xaxt = 'n', yaxt = 'n')
    }
    for (i in 1:nrow(map)) {
      points(x = map[i, 1], y = map[i, 2], pch = 20, col = hcol[i], cex = cexbg[i])
      points(x = map[i, 1], y = map[i, 2], pch = 20, col = cols[i], cex = cex[i])
      if (!is.null(show_cluster_number)) {
        text(x = map[i, 1], y = map[i, 2], labels = srt[i], adj = 0.5, cex = show_cluster_number)
      }
    }
  } else {
    plot(map, bty = bty, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20, col = cols, cex = cex, main = tit)
    if (contours) {
      kde = bkde2D(panoply_obj$map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
      contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE, col = 'darkgray', add = TRUE, xaxt = 'n', yaxt = 'n')
    }
  }
}

plot_tsne_spread = function(panoply_obj, markers = NULL, mode = c("arithmetic", "robust"), cex = 50.0, proportional = TRUE, emph = TRUE, highlight_clusters, contours = FALSE, show_cluster_number = NULL) {
  mode = match.arg(mode)

  # calculate plot layout
  n = length(markers) + 1
  sq = sqrt(n)
  frac = sq - floor(sq)
  if (frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  par(mfrow = c(dn, ac), mar = c(1, 1, 2, 1))
  for (i in 1:length(markers)) {
    plot_tsne(panoply_obj = panoply_obj, marker = markers[i], mode = mode, cex = cex, proportional = proportional, emph = emph, highlight_clusters = highlight_clusters, contours = contours, show_cluster_number = show_cluster_number)
  }
}


