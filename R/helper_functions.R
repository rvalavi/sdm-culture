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
    
    warning("Number of duplicate records removed: ", sum(c(!pr_dup, !bg_dup)), "\n")
    
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
    
    return(
        if (drop_na) {
            tidyr::drop_na(vals)
        } else {
            vals
        }
    )
}
