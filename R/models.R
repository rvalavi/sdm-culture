# install required packages
pkg <- c(
    "ranger",
    "dismo",
    "gbm",
    "mgcv",
    "glmnet",
    "terra",
    "blockCV",
    "dismo",
    "precrec",
    "caret"
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
    library(dismo) #***** or predicts
    library(precrec)
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
ensemble <- function(
        x,
        y = "",
        fold_ids = NULL,
        models = c("GLM", "GAM", "GBM", "RF")
) {
    
    # get covariate names
    covars <- names(x)
    covars <- covars[covars != y]
    
    if(is.null(fold_ids)){
        fold_ids <- dismo::kfold(x, k = 5)
    }
    
    # if null use all
    if (is.null(models)) models <- c("GLM", "GAM", "GBM", "RF")
    if (tolower(models[1]) == "all") models <- c("GLM", "GAM", "GBM", "RF")
    # set the default to null
    ls_mod <- gm_mod <- br_mod <- rf_mod <- NULL
    
    # compute weights
    prNum <- as.numeric(table(x[, y, drop=TRUE])["1"]) # number of presences
    bgNum <- as.numeric(table(x[, y, drop=TRUE])["0"]) # number of backgrounds
    wt <- ifelse(x[, y, drop = TRUE] == 1, 1, prNum / bgNum)
    
    if ("GAM" %in% models) {
        require(mgcv)
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
            family = binomial(link = "logit"),
            weights = wt,
            method = "REML",
            select = TRUE
        )
    }
    
    if ("RF" %in% models) {
        message("Fitting RF...")
        
        warning("Fitting RF needs to consider down-sampling!")
        
        # fit random forest with ranger
        rf_mod <- random_forest(
            data = x,
            y = y,
            max.depth = c(0, 5, 10),
            splitrule = c("gini", "hellinger", "extratrees"),
            num.trees = c(1000),
            mtry = 2:4,
            foldid = fold_ids,
            threads = NULL,
            plot = TRUE
        )
    }
    
    if ("GBM" %in% models) {
        message("Fitting GBM...")
        # fit boosted regression tress (BRT) a.k.a GBM
        br_mod <- dismo::gbm.step(
            data = x,
            gbm.y = which(names(x) == y),
            gbm.x = which(names(x) != y),
            family = "bernoulli",
            tree.complexity = ifelse(prNum < 50, 1, 5),
            learning.rate = 0.001,
            bag.fraction = 0.75,
            max.trees = 10000,
            n.trees = 50,
            site.weights = wt,
            fold.vector = fold_ids,
            n.folds = length(unique(fold_ids)),
            silent = TRUE
        )
    }
    
    if ("GLM" %in% models) {
        message("Fitting GLM-Lasso...")
        # fit regularized regression with L1
        quad_obj <- make_quadratic(x, cols = covars)
        training_quad <- predict(quad_obj, newdata = x)
        new_vars <- names(training_quad)[names(training_quad) != y]
        training_sparse <- sparse.model.matrix(~., training_quad[, new_vars])
        
        ls_mod <- glmnet::cv.glmnet(
            x = training_sparse,
            y = training_quad[, y],
            family = "binomial",
            alpha = 1, # fitting lasso
            weights = wt,
            foldid = fold_ids
        )
        plot(ls_mod)
    }
    
    out <- list(
        "GLM-Lasso" = ls_mod,
        "GAM" = gm_mod,
        "GBM" = br_mod,
        "RF" = rf_mod,
        "quad" = quad_obj
    )
    class(out) <- c("ensemble", "list")
    
    return(out)
}

# print method for ensemble object
print.ensemble <- function(x, ...){
    print(class(x))
    mods <- names(x)[!sapply(x, is.null)]
    mods <- mods[mods != "quad"]
    cat("Models:", paste(mods, collapse = ", "))
}

# predict with ensemble object
predict.ensemble <- function(object, newdata, ...){
    require(ranger)
    require(dismo)
    require(gbm)
    require(mgcv)
    require(glmnet)
    
    pred_ls <- pred_gm <- pred_rf <- pred_br <- NA
    
    # predict lasso with the 1 sd lambda
    if (!is.null(object[["GLM-Lasso"]])) {
        pred_quad <- predict(object[["quad"]], newdata = newdata)
        pred_sparse <- sparse.model.matrix( ~., pred_quad)
        pred_ls <- predict(object[["GLM-Lasso"]], pred_sparse, s = "lambda.min")
    }
    
    # predict gam model
    if (!is.null(object[["GAM"]])) {
        pred_gm <- predict(object[["GAM"]], newdata)
    }
    
    # predict rf with ranger
    if (!is.null(object[["RF"]])) {
        pred_rf <- predict(object[["RF"]], newdata)$predictions
    }
    
    # predict brt with the best tree number
    brt <- object[["GBM"]]
    if (!is.null(brt)) {
        pred_br <- predict(brt, newdata, n.trees = brt$gbm.call$best.trees)
    }
    
    out <- data.frame(
        l = pred_ls[, 1],
        g = pred_gm,
        r = pred_rf,
        b = pred_br
    )
    
    return(
        rowMeans(pred_ls, na.rm = TRUE)
    )
}

# function for simultaneous tuning of maxent regularization multiplier and features
maxent_param <- function(data, y = "occ", k = 5, folds = NULL, filepath){
    require(dismo)
    require(caret)
    require(precrec)
    if(is.null(folds)){
        # generate balanced CV folds
        folds <- caret::createFolds(y = as.factor(data$occ), k = k)
    }
    names(data)[which(names(data) == y)] <- "occ"
    covars <- names(data)[which(names(data) != y)]
    # regularization multipliers
    ms <- c(0.5, 1, 2, 3, 4)
    grid <- expand.grid(
        regmult = paste0("betamultiplier=", ms),
        features = list(
            c("noautofeature", "nothreshold"), # LQHP
            c("noautofeature", "nothreshold", "noproduct"), # LQH
            c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
            c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
            c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
        stringsAsFactors = FALSE
    )
    AUCs <- c()
    for(n in seq_along(grid[,1])){
        full_pred <- data.frame()
        for(i in seq_len(length(folds))){
            trainSet <- unlist(folds[-i])
            testSet <- unlist(folds[i])
            if(inherits(try(
                maxmod <- dismo::maxent(x = data[trainSet, covars],
                                        p = data$occ[trainSet],
                                        removeDuplicates = FALSE,
                                        path = filepath,
                                        args = as.character(unlist(grid[n, ]))
                )
            ), "try-error")){
                next
            }
            modpred <- predict(maxmod, data[testSet, covars], args = "outputformat=cloglog")
            pred_df <- data.frame(score = modpred, label = data$occ[testSet])
            full_pred <- rbind(full_pred, pred_df)
        }
        AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score, 
                                                 labels = full_pred$label))[1,4]
    }
    best_param <- as.character(unlist(grid[which.max(AUCs), ]))
    return(best_param)
}


# tune parameters in ranger mode
random_forest <- function(data, 
                          y = "occ", 
                          max.depth = c(0, 5, 10, 15),
                          splitrule = c("gini", "hellinger"),
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
    
    return(
        outCells[!is.na(outCells)]
    )
}


# Orthogonal quadratic polynomials for glmnet
#
# A function to create quadratic terms for glmnet functions i.e. lasso and ridge regression.
# The output is an object of make_quadratic that can be used to predict on rasters and data.frames
# for creating the quadratic terms.
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

