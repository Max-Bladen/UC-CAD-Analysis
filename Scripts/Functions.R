
### ======================================================================== ###
### Load libraries                                                           ###
### ======================================================================== ###

suppressMessages(library(gtools))
suppressMessages(library(mixOmics))
suppressMessages(library(readxl))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(UniprotR))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyr))
suppressMessages(library(BiocParallel))


### ======================================================================== ###
### Global variables                                                         ###
### ======================================================================== ###

DEFAULT.PLOT.WIDTH.MM=150
DEFAULT.PLOT.HEIGHT.MM=100





### ======================================================================== ###
### Data manipulation functions                                              ###
### ======================================================================== ###

### ======================================================================== ###
#
#
Sort_Loadings <- function(obj,
                          block = NULL,
                          comp=1) {
  if (is.null(block)) { loadings <- obj$loadings[,comp] }
  else { loadings <- obj$loadings[[block]][,comp] }
  
  order <- sort(abs(loadings), decreasing=T, index.return=T)$ix
  return(loadings[order])
}



### ======================================================================== ###
# Takes the sample name and extracts which group it belongs to
#
Sample_to_Group <- function(sample) {
  # Iterate over all items in input vector/list
  unlist(lapply(sample, function(s) {
    if (length(grep("noCAD", s))==1) { return("noCAD") }
    else if (length(grep("CAD", s))==1) { return("CAD") }
    # Return NA if cannot find either noCAD or CAD
    return(NA)
  }))
}



### ======================================================================== ###
# Takes a lipid feature name and extracts which lipid class it belongs to
#
Lipid_to_Class <- function(molecule) {
  # Iterate over all items in input vector/list
  unlist(lapply(molecule, function(m) {
    # If it starts with TAG, return immediately
    if (substr(m, 1, 3) == "TAG") {return("TAG")}
    # Otherwise, remove redundant characters and return remainder
    gsub("\\(.*$", "", m)
  }))
  
}



### ======================================================================== ###
# Converts a data.table column into a vector
#
v <- function(dtCol) {
  colClass <- class(as.vector(as.matrix(dtCol)))
  v <- as.vector(as.matrix(dtCol))
  class(v) <- colClass
  v
}



### ======================================================================== ###
# Cacluates the absolute correlation matrix of the columns of a given data table
#
AbsCorMat <- function(df) {
  # Calculate correlations between columns, find absolute value and conert
  # to data.table
  corMat <- df %>% cor() %>% abs() %>% as.data.table()
  # Clean up unncessary values
  corMat[lower.tri(corMat)] <- NA
  diag(corMat) <- NA
  # return
  corMat
}



### ======================================================================== ###
# Takes a dataframe and a maximum intercorrelation threshold. Iteratively 
# determine which feature is contained in the highest intercorrelation pair,
# remove it and repeat. Continues until there are no pairs with a correlation
# above threshold
# 
Remove_Intercorrelated_Features <- function(df,
                                            threshold = 0.9) {
  # Initial absolute correlation matrix 
  corMat <- AbsCorMat(df)
  # What is maximum value
  maxCor <- corMat %>% max(na.rm=T)
  
  # Iterate
  while (maxCor > threshold) {
    # Determine which feature has this highest value
    highCorFeat <- which(corMat == max(corMat, na.rm=T), 
                         arr.ind=T)[2]
    
    # Remove it from dataframe and corMat
    df <- df[,-highCorFeat, with = F]
    corMat <- corMat[-highCorFeat,-highCorFeat, with = F]
    
    # Re-determine maximum value
    maxCor <- corMat%>% max(na.rm=T)
  }
  
  df
}



### ======================================================================== ###
# Iterates over a range of intercorrelation thresholds, calculates the number
# of retained features (using the Remove_Intercorrelated_Features function)
# and then plots these values as a line plot.
#
Plot_Intercorrelation_Threshold <- function(df, 
                                            range = seq(0.1, 1, 0.05),
                                            title = NULL,
                                            showPlot = F,
                                            savePath = NULL) {
  # For each value in the range parameter, remove intercorrelated features
  # Then count the number of features in each of these
  feats <- lapply(range, function(t){
    Remove_Intercorrelated_Features(df, t)
  }) %>% lapply(ncol) %>% unlist()
  
  # Plot
  plot <- data.table(threshold = range,
             feats = feats) %>% 
    ggplot(aes(x = threshold,
               y = feats)) + 
    scale_x_reverse() + 
    geom_line() + geom_point() + 
    xlab("Intercorrelation threshold") + 
    ylab("Number of retained features") + 
    theme_bw()
  
  # Add title if desired
  if (!is.null(title)) {
    plot <- plot + ggtitle(title)
  }
  
  Save_And_Plot(plot, showPlot, savePath)
}



### ======================================================================== ###
# 
#
Matrix_Var <- function(x, dim = 1) {
  if(dim == 1){
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  } else if (dim == 2) {
    rowSums((t(x) - colMeans(x))^2)/(dim(x)[1] - 1)
  } else stop("Please enter valid dimension")
}





### ======================================================================== ###
### Plotting functions                                                       ###
### ======================================================================== ###

### ======================================================================== ###
# General, all-purpose colour palette
Palette <- colorRampPalette(c("#881199", "#881111","#ED8C32","#117733","#6699CC"))

# Colours for sample groups
Group.Palette <- c("#881199", "#6699CC")

# Colours for plots using two blocks
Block.Palette <- c("#2c03fc", "#315d00")

# Colours for plots using three blocks
Block.Palette.Extra <- c("#2c03fc", "#315d00", "#7f0b9c")



### ======================================================================== ###
# Generalised function to save plot to file and return plot if desired.
# 
Save_And_Plot <- function(plot = NULL,
                          showPlot = F,
                          savePath = NULL,
                          make.dir = F,
                          w = DEFAULT.PLOT.WIDTH.MM, 
                          h = DEFAULT.PLOT.HEIGHT.MM,
                          units = "mm",
                          verbose = F) {
  # show the plot if desired
  if (showPlot) {print(plot)}
  
  # save plot if desired
  if (!is.null(savePath)) {
    if ((!dir.exists(dirname(savePath)) & make.dir) | dir.exists(dirname(savePath))) {
      if (make.dir) {message("Making directory...")}  
      ggsave(savePath, plot,
             width = w,
             height = h,
             units = units) 
    }
    else {message("Directory not found!")}
    if (verbose) {message("Plot saved to :"); message(paste0("    ", savePath))}
  }
  
  plot
} 



### ======================================================================== ###
# Takes a list of ggplot objects and combines them into a grid plot. Various
# options for adjusting formatting and style
#
Plot_Grid <- function(plot.list = NULL,
                      labels = NULL,
                      n.row = NULL,
                      n.col = NULL,
                      widths = NULL,
                      heights = NULL,
                      title = NULL,
                      shared.x.lab = NULL,
                      shared.y.lab = NULL,
                      shared.lgnd.pos = NULL,
                      w = NULL, h = NULL,
                      w.pp = DEFAULT.PLOT.WIDTH.MM,
                      h.pp = DEFAULT.PLOT.HEIGHT.MM,
                      showPlot = F,
                      savePath = NULL) {
  
  Factors <- function(x,
                      prime.factors=T) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    if (!prime.factors & length(factors)>2) {
      factors <- factors[factors!=1]
      factors <- factors[factors!=x]
    }
    return(factors)
  }
  
  Determine.Plot_Grid.Dims <- function(N) {
    sqrt.N <- sqrt(N)
    if (floor(sqrt.N) == sqrt.N) {nrow <- ncol <- sqrt.N } 
    else {
      factors <- Factors(N, F)
      diff_factors <- diff(sort(as.numeric(factors)))
      min_diff <- min(diff_factors)
      closest_factors <- sort(as.numeric(factors))[which(diff_factors == min_diff)]
      
      if (length(closest_factors) == 1) {
        nrow <- ceiling(N / closest_factors)
        ncol <- closest_factors
      } else {
        nrow <- closest_factors[1]
        ncol <- closest_factors[2]
      }
    }
    return(c(nrow, ncol))
  }
  
  n.p <- length(plot.list)
  if (is.null(n.row) & is.null(n.col)) {
    dims <- Determine.Plot_Grid.Dims(n.p)
    n.row <- dims[1]
    n.col <- dims[2]
    rm(dims)
  }
  if (is.null(w)) {w <- n.row*w.pp}
  if (is.null(h)) {h <- n.row*h.pp}
  
  if (!is.null(labels)) {
    if (is.logical(labels)) {
      if (labels) {
        labels <- paste0(LETTERS, ")")[1:n.p]
      } else {labels <- NULL}
    }
  }
  
  
  
  if (is.null(widths)) {widths <- rep(1, n.p)}
  if (is.null(heights)) {heights <- rep(1, n.p)}
  rm(n.p)
  
  if (!is.null(shared.x.lab)) {
    plot.list <- lapply(plot.list, function(p) {p + rremove("xlab")})
  }
  if (!is.null(shared.y.lab)) {
    plot.list <- lapply(plot.list, function(p) {p + rremove("ylab")})
  }
  
  if (!is.null(shared.lgnd.pos)) {
    combined <- ggarrange(plotlist = plot.list,
                          widths = widths,
                          heights = heights,
                          labels = labels, 
                          common.legend = T, legend = shared.lgnd.pos,
                          nrow=n.row, ncol=n.col)
  }
  else {
    combined <- ggarrange(plotlist = plot.list,
                          widths = widths,
                          heights = heights,
                          labels = labels, 
                          common.legend = F, 
                          nrow=n.row, ncol=n.col)
  }
  
  if (!is.null(title)) {
    combined <- annotate_figure(combined,
                                fig.lab = title,
                                fig.lab.pos = "top.left",
                                fig.lab.size = 15,
                                top="")
  }
  
  if (!is.null(shared.x.lab)) {
    combined <- annotate_figure(combined,
                                bottom = text_grob(shared.x.lab))
  }
  if (!is.null(shared.y.lab)) {
    combined <- annotate_figure(combined,
                                left = text_grob(shared.y.lab, 
                                                 rot = 90))
  }
  combined <- combined + bgcolor("white")
  
  Save_And_Plot(combined, showPlot, savePath, w=w, h=h)
}



### ======================================================================== ###
### Misc functions                                                           ###
### ======================================================================== ###

### ======================================================================== ###
# Opposite of '%in%' operator
'%!in%' <- function(x,y)!('%in%'(x,y))



### ======================================================================== ###
# Equivalent of the base `head()` function but controls the number of columns
# shown. This is purely a QoL exploration function and was not used in analysis
#
Head <- function(df,
                 nr = 10,
                 nc = 5) {
  # Get the number of rows and columns in the input data frame
  r <- nrow(df)
  c <- ncol(df)
  
  # Ensure that nc and nr are within the bounds of the data frame dimensions
  if (nc > c) {
    nc <- c
  }
  if (nr > r) {
    nr <- r
  }
  
  # Extract the top-left portion of the data frame
  return(df[1:nr, 1:nc])
}



### ======================================================================== ###
# Equivalent to the base::table() function, but includes any values with a count
# of 0. 
# 
Full.Table <- function(vec,
                       all.vals=unique(vec)) {
  vals <- unique(vec)
  all.vals <- as.character(all.vals)
  # use table() to find frequency of items in vector
  tbl <- table(unname(vec))
  # for any missing vals from all.vals, add them to tbl with a count of 0
  for (val in all.vals[all.vals %!in% vals]) { tbl[val] <- 0 }
  # order properly
  tbl <- tbl[mixedorder(names(tbl))]
  # remove any vals not in all.vals
  tbl <- tbl[all.vals]
  
  return(as.table(tbl))
}



### ======================================================================== ###
# For a given confusion matrix, calculate the sensitivity, specificity, PPV and
# NPV
#
Extract_Confusion_Metrics <- function(confMat) {
  
  TP <- confMat[1,1]
  TN <- confMat[2,2]
  FP <- confMat[2,1]
  FN <- confMat[1,2]
  
  sens <- (TP / (TP + FN)) * 100
  spec <- (TN / (TN + FP)) * 100
  ppv <- (TP / (TP + FP)) * 100
  npv <- (TN / (TN + FN)) * 100
  
  return(list(Sensitivity=sens,
              Specificity=spec,
              PPV=ppv,
              NPV=npv))
}



### ======================================================================== ###
# Used to separate a data object into two (similarly structure) data objects
# by the provided sample indices
# 
Extract_Data_Subset <- function(data,
                                idx) {
  # Initialise returned objects
  newData <- oldData <- data
  
  # Set up new data object to contain provided samples
  newData$samples <- newData$samples[idx,]
  newData$dfs <- newData$dfs %>% lapply(function(df) {
    df[idx,]
  })
  
  # Set up old data object to contain all samples except those provided
  oldData$samples <- oldData$samples[-idx,]
  oldData$dfs <- oldData$dfs %>% lapply(function(df) {
    df[-idx,]
  })
  
  return(list(oldData = oldData,
              newData = newData))
  
}