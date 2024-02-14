# install required packages
pkg <- c(
  "ranger",
  "dismo",
  "gbm",
  "mgcv",
  "glmnet",
  "terra",
  "blockCV",
  "MASS",
  "Rcpp"
)
pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
if(length(pkgna)){
  nm <- paste(pkgna, collapse = ", ")
  message("Installing packages: ", nm, "...")
  utils::install.packages(pkgna)
}

# remotes::install_github("rvalavi/myspatial")
# library(myspatial)
suppressPackageStartupMessages({
  library(ranger)
  library(dismo)
  library(gbm)
  library(mgcv)
  library(glmnet)
})


#' Fitting several models as an ensemble
#'
#' A function to ... 
#' The response if Gaussin here..
#' Models include: GAM, Lasso, RF and BRT (a.k.a. GBM)
#'
#' @param x data.frame of the training data.
#' @param y character, the column name of the response
#'
#' @author Roozbeh Valavi
#'
#' @return an object of ensemble containing individual model objects for prediction
#' @export
#'
#' @examples
ensemble <- function(x, y = "", fold_ids = NULL){
  # get covariate names
  covars <- names(x)
  covars <- covars[covars != y]
  
  if(is.null(fold_ids)){
    fold_ids <- dismo::kfold(x, k = 5)
  }
  
  message("Fitting GAM...")
  # fit generalised additive model (GAM))
  gm_form <- paste(
    y, 
    paste0("s(", covars, ")", collapse = " + "),
    sep = " ~ "
  )
  gm_mod <- mgcv::gam(
    formula = as.formula(gm_form),
    data = x,
    select = TRUE,
    method = "REML"
  )
  
  message("Fitting RF...")
  # fit random forest with ranger
  rf_mod <- random_forest(
    data = x,
    y = y,
    max.depth = c(0, 5, 10),
    splitrule = c("variance", "extratrees", "maxstat"),
    num.trees = c(500),
    mtry = 2:4,
    foldid = fold_ids,
    threads = NULL,
    plot = TRUE
  )

  message("Fitting BRT...")
  # fit boosted regression tress (BRT) a.k.a GBM
  br_mod <- dismo::gbm.step(
    data = x,
    gbm.y = which(names(x) == y),
    gbm.x = which(names(x) != y),
    family = "gaussian",
    tree.complexity = 5,
    learning.rate = 0.005,
    bag.fraction = 0.75,
    max.trees = 10000,
    n.trees = 50,
    fold.vector = fold_ids,
    n.folds = length(unique(fold_ids)),
    silent = TRUE
  )
  
  message("Fitting Lasso...")
  # fit regularised regression with L1
  quad_obj <- make_quadratic(x, cols = covars)
  training_quad <- predict(quad_obj, newdata = x)
  new_vars <- names(training_quad)[names(training_quad) != y]
  training_sparse <- sparse.model.matrix(~., training_quad[, new_vars])
  ls_mod <- glmnet::cv.glmnet(
    x = training_sparse,
    y = training_quad[, y],
    family = "gaussian",
    alpha = 1, # fitting lasso
    foldid = fold_ids
    # nfolds = 5
  )
  plot(ls_mod)
  
  out <- list(
    "gam" = gm_mod,
    # "rf" = rf_mod,
    # "brt" = br_mod,
    "lasso" = ls_mod,
    "quad" = quad_obj
  )
  class(out) <- c("ensemble", "list")
  
  return(out)
}

#' @export
#' @method print ensemble
print.ensemble <- function(x, ...){
  print(class(x))
}

# predict with ensemble object
#' @export
#' @method predict ensemble
predict.ensemble <- function(object, newdata, ...){
  require(ranger)
  require(dismo)
  require(gbm)
  require(mgcv)
  require(glmnet)
  # predict gam model
  pred_gm <- predict(object[["gam"]], newdata)
  # predict rf with ranger
  pred_rf <- predict(object[["rf"]], newdata)$predictions
  # predict brt with the best tree number
  brt <- object[["brt"]]
  pred_br <- predict(brt, newdata, n.trees = brt$gbm.call$best.trees)
  # predict lasso with the 1 sd lambda
  pred_quad <- predict(object[["quad"]], newdata = newdata)
  pred_sparse <- sparse.model.matrix( ~., pred_quad)
  pred_ls <- predict(object[["lasso"]], pred_sparse, s = "lambda.min")
  
  out <- data.frame(
    g = pred_gm,
    r = pred_rf,
    b = pred_br,
    l = pred_ls[, 1]
  )
  
  return(
    rowMeans(pred_ls)
  )
}


# tune parameters in ranger mode
random_forest <- function(data, 
                          y = "occ", 
                          max.depth = c(0, 5, 10, 15),
                          splitrule = c("variance", "extratrees"),
                          num.trees = 1000,
                          mtry = NULL,
                          foldid = NULL,
                          threads = NULL,
                          plot = TRUE){
  
  require(ranger)
  require(ggplot2)
  
  form <- as.formula(
    paste(y, "~ .")
  )
  
  if(is.null(mtry)){
    mtry <- floor(sqrt(ncol(data) - 1))
  }
  
  grid <- expand.grid(max.depth = max.depth, 
                      splitrule = splitrule,
                      num.trees = num.trees,
                      mtry = mtry,
                      stringsAsFactors = FALSE)
  
  if(is.null(foldid)){
    foldid <- dismo::kfold(data, k = 5)
  }
  nfold <- sort(unique(foldid))
  
  evalmodel <- data.frame(depth = rep(NA, nrow(grid)), split = NA)
  for(i in seq_along(grid[,1])){
    modrsq <- c()
    for(k in nfold){
      
      train_set <- which(foldid != k)
      test_set <- which(foldid == k)
      
      mod <- ranger::ranger(
        formula = form,
        data = data[train_set, ], 
        num.trees = grid$num.trees[i],
        splitrule = grid$splitrule[i],
        max.depth = grid$max.depth[i],
        mtry = grid$mtry[i],
        num.threads = threads,
        replace = TRUE
      )
      
      pred <- predict(mod, data[test_set, ])$predictions
      
      modrsq[k] <- cor(pred, data[test_set, y, drop = TRUE])^2
    }
    evalmodel$depth[i] <- grid$max.depth[i]
    evalmodel$split[i] <- grid$splitrule[i]
    evalmodel$ntrees[i] <- grid$num.trees[i]
    evalmodel$mtry[i] <- grid$mtry[i]
    evalmodel$eval[i] <- mean(modrsq)
    # evalmodel$aucse[i] <- sd(modrsq) / sqrt(nfold)
  }
  
  bestparam <- which.max(evalmodel$eval)
  
  if(plot){
    print(
      plot_tune_param(evalmodel)
    )
  }
  print(evalmodel[bestparam, ])
  
  finalmod <- ranger::ranger(
    formula = form,
    data = data, 
    num.trees = evalmodel$ntrees[bestparam],
    mtry = evalmodel$mtry[bestparam],
    splitrule = evalmodel$split[bestparam],
    max.depth = evalmodel$depth[bestparam],
    num.threads = threads,
    replace = TRUE
  )
  
  return(finalmod)
}

# plot tuning parameters of ranger
plot_tune_param <- function(x){
  levs <- sort(unique(x$depth))
  x$depth <- as.character(x$depth)
  x$depth <- ifelse(x$depth == 0, "max", x$depth)
  x$depth <- factor(x$depth, levels = c(as.character(levs[-1]), "max"))
  x$id <- as.numeric(as.factor(x$mtry))
  return(
    ggplot(data = x, aes(x = id, y = eval, col = depth, shape = split)) +
      geom_point(aes(size = as.factor(mtry))) +
      geom_path() +
      facet_grid(~ntrees) +
      scale_x_continuous(breaks = unique(x$id)) +
      theme_bw() +
      labs(x = "", y = "R-squared", size = "mtry")
  )
}


# function to remove the training cells and their direct neighbours
cell_neighbour <- function(r, cell, size = 3){
  h <- (size - 1) / 2
  dif <- rep(seq(-h, h), each = size)
  rowcols <- rowColFromCell(r, cell = cell)
  newr <- matrix(dif, nrow = size, byrow = TRUE) + rowcols[1]
  newc <- matrix(dif, nrow = size, byrow = FALSE) + rowcols[2]
  outCells <- cellFromRowCol(r, row = as.vector(newr), col = as.vector(newc))
  return(outCells[!is.na(outCells)])
}

# clean climate variable name
simplify_name <- function(x){
  str_split(x, "_") %>% 
    map(pluck, 2) %>% 
    unlist()
}


#' Orthogonal quadratic polynomials for glmnet
#'
#' A function to create quadratic terms for glmnet functions i.e. lasso and ridge regression.
#' The output is an object of make_quadratic that can be used to predict on rasters and data.frames
#' for creating the quadratic terms.
#'
#' @param df a data.frame, typically the training data.
#' @param cols the name or index of the columns to be transformed. If NULL, all the columns will be transformed.
#' The factor columns won't be transformed.
#' @param verbose logical. print messages.
#'
#' @author Roozbeh Valavi
#'
#' @return an object of make_quadratic that can be used to predict on rasters and data.frames
#' @export
#'
#' @examples
make_quadratic <- function(df, cols = NULL, verbose = TRUE){
  if(is.null(cols)){
    cols <- colnames(df)
  }
  if(is.numeric(cols)){
    cols <- colnames(df)[cols]
  }
  # remove the factors
  if(any(sapply(df[,cols], is.factor))){
    if(verbose){
      message("The factor columns were removed form cols: ", cols[which(sapply(df[,cols], is.factor))])
    }
    cols <- cols[-which(sapply(df[,cols], is.factor))]
  }
  if(!all(is.element(cols, colnames(df)))){
    stop("The cols should be the same as the column names.")
  }
  xbar <- apply(df[,cols], 2, mean)
  x1 <- data.frame(mapply(`-`, df[,cols], xbar, SIMPLIFY = FALSE))
  alpha <- colSums(x1 ^ 3) / colSums(x1 ^ 2)
  # specify the output class
  finalList <- list(names = cols, xbars = xbar, alphas = alpha)
  class(finalList) <- c("make_quadratic")
  return(finalList)
}

#' @export
#' @method predict make_quadratic
predict.make_quadratic <- function(object, newdata, ...){
  if(!methods::is(object, "make_quadratic"))
    stop("object should be a make_quadratic object.")
  if(!all(object$names %in% names(newdata)))
    stop("The newdata does not have the same names as the object.")
  ncl <- object$names
  if(methods::is(newdata, "Raster")){
    for(i in ncl){
      x1 <- newdata[[i]] - object$xbars[i]
      x2 <- (x1 ^ 2) - (object$alphas[i] * x1)
      if(raster::nlayers(newdata) > 1){
        newdata <- newdata[[-which(names(newdata) == i)]]
        newdata <- raster::stack(newdata, x1)
      } else{
        newdata <- x1
      }
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_1")
      newdata <- raster::stack(newdata, x2)
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_2")
    }
  } else if(methods::is(newdata, "data.frame")){
    for(i in ncl){
      x1 <- newdata[,i] - object$xbars[i]
      x2 <- x1 ^ 2 - object$alphas[i] * x1
      newdata <- newdata[,-which(names(newdata) == i), drop = FALSE]
      newdata[,ncol(newdata) + 1] <- x1
      names(newdata)[ncol(newdata)] <- paste0(i, "_1")
      newdata[,ncol(newdata) + 1] <- x2
      names(newdata)[ncol(newdata)] <- paste0(i, "_2")
    }
  } else stop("newdata should be a raster or a data.frame.")
  return(newdata)
}
