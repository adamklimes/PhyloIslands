# adapted from "diversitree::trait.plot"
trait.plot2 <- function (tree, dat, cols, lab = names(cols), str = NULL, class = NULL, 
    type = "f", w = 1/50, legend = length(cols) > 1, cex.lab = 0.5, 
    font.lab = 3, cex.legend = 0.75, margin = 1/4, check = TRUE, 
    quiet = FALSE, stl = TRUE, g_lwd = 1.5, ...) 
{
    if (!(type %in% c("f", "p"))) 
        stop("Only types 'f'an and 'p'hylogram are available")
    if (!is.null(class) && length(class) != length(tree$tip.label)) 
        stop("'class' must be a vector along tree$tip.label")
    n <- length(cols)
    if (n < 1) 
        stop("Need some colours")
    if (!is.data.frame(dat)) {
        if (is.vector(dat) && n == 1) {
            nm <- names(dat)
            dat <- matrix(dat)
            rownames(dat) <- nm
        }
        else {
            stop("dat must be a matrix")
        }
    }
    if (!all(tree$tip.label %in% rownames(dat))) 
        stop("All taxa must have entries in 'dat' (rownames)")
    if (n > 1) {
        if (!all(names(cols) %in% names(dat))) 
            stop("Not all colours have data")
        if (is.null(names(cols))) 
            stop("'cols' must be named")
        dat <- dat[names(cols)]
    }
    if (is.null(str)) {
        str <- lapply(dat, function(x) as.character(sort(unique(x))))
    }
    dat <- dat[tree$tip.label, , drop = FALSE]
#    par(mar = rep(0, 4))
    t <- max(branching.times(tree))
    w <- w * t
    if (is.null(class)) {
        plt <- diversitree:::plot2.phylo(tree, type = type, show.tip.label = stl, 
            label.offset = (n + 2) * w, cex = cex.lab, ...)
    }
    else {
        plt <- diversitree:::plot2.phylo(tree, type = type, show.tip.label = FALSE, 
            label.offset = t * margin, ...)
        diversitree:::group.label.tip(plt, class, c("grey45", "black"), "black", 
            lwd = g_lwd, offset.bar = w * (n + 3), offset.lab = w * 
                (n + 4), cex = cex.lab, font = font.lab, check = check, 
            quiet = quiet)
    }
    if (type == "f") {
        xy <- plt$xy
        theta <- xy$theta[seq_along(tree$tip.label)]
        dt <- diff(sort(theta))[1]/2
        for (i in seq_along(cols)) {
            idx <- dat[[names(dat)[i]]]
            if (any(idx == 0, na.rm = TRUE)) 
                idx <- idx + 1
            diversitree:::filled.arcs(theta - dt, theta + dt, max(xy$x) + i * 
                w, w, cols[[i]][idx])
        }
    }
    else {
        xy <- plt$xy[seq_along(tree$tip.label), ]
        dy <- 0.5
        for (i in seq_along(cols)) {
            idx <- dat[[names(dat)[i]]]
            if (any(idx == 0, na.rm = TRUE)) 
                idx <- idx + 1
            xleft <- xy[1, 1] + w * i
            xright <- xleft + w
            ybottom <- xy[, 2] - dy
            ytop <- ybottom + dy * 2
            rect(xleft, ybottom, xright, ytop, col = cols[[i]][idx], 
                border = NA)
        }
    }
    if (legend) {
        for (i in seq_along(cols)) {
            c.i <- cols[[i]]
            leg.txt <- str[[i]]
            leg.arg <- list(legend = leg.txt, title = lab[i], 
                title.adj = 0, bty = "n", fill = c.i, cex = cex.legend, 
                horiz = TRUE)
            ifelse(i == 1, leg <- do.call("legend", c("topleft", 
                leg.arg)), leg <- do.call("legend", c(leg$rect$left, 
                leg$rect$top - leg$rect$h, leg.arg)))
        }
    }
    invisible(plt)
}