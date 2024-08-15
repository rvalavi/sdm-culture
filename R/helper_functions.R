# filter coordinates based on cells
# NOTE: this function will change the order of the rows
rm_duplicates <- function(x, r, column = "occ") {
    # filter presences and background samples to remove duplicates separately
    pr <- dplyr::filter(x, !!rlang::sym(column) == 1)
    bg <- dplyr::filter(x, !!rlang::sym(column) == 0)
    # find the duplicate points in a pixel
    pr_dup <- duplicated(
        terra::cellFromXY(r[[1]], as.matrix(pr[, c("x", "y")]))
    )
    bg_dup <- duplicated(
        terra::cellFromXY(r[[1]], as.matrix(bg[, c("x", "y")]))
    )
    
    num_rm <- sum(c(!pr_dup, !bg_dup))
    if (num_rm > 0)
        warning("Number of duplicate records removed: ", num_rm, "\n")
    
    return(
        rbind(
            pr[!pr_dup, ],
            bg[!bg_dup, ]
        )
    )
}

# extract values and drop na
extract_value <- function(r, x, drop_na = TRUE) {
    vals <- cbind(
        x,
        terra::extract(r, x[, c("x", "y")], ID = FALSE)
    )
    
    ccs <- complete.cases(vals)
    num_na <- sum(!ccs)
    
    return(
        if (drop_na) {
            warning("Number of NA rows removed: ", num_na, "\n")
            vals[ccs, ]
        } else {
            if (num_na > 0) warning("There ", num_na, " NAs rows.\n")
            vals
        }
    )
}
