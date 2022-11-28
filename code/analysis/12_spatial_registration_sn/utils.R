## Polish annotations
fix_layer_order <- function(l) {
    star <- ifelse(grepl("\\*", l), "*", "")
    l <- gsub("L|\\*", "", l)
    l <- sort(unlist(strsplit(l, "/")))

    if (all(l == "WM")) {
        return(paste0("WM", star))
    }

    l[[1]] <- paste0("L", l[[1]])
    # if ("WM" %in% l) l <- c("WM", l[l != "WM"])
    fix <- paste0(paste0(l, collapse = "/"), star)

    return(fix)
}

fix_layer_order2 <- Vectorize(fix_layer_order)
