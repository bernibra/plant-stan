library('igraph')
library("TeachingDemos")
library("grid")
library("extrafont")
library(cluster, lib.loc = .Library)
library("factoextra")
library("magrittr")
library("NbClust")
library(pracma)
library(vegan)
library(ggplot2)
library(amap)


`plot.betadisper` <- function(x, axes = c(1,2), cex = 0.7, pch = seq_len(ng),
                              col = NULL, lty = "solid", lwd = 1, hull = TRUE,
                              ellipse = FALSE, conf,
                              segments = TRUE, seg.col = "grey",
                              seg.lty = lty, seg.lwd = lwd,
                              label = TRUE, label.cex = 1,
                              ylab, xlab, main, sub, ...)
{
  localAxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
  localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
  localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
  Ellipse <- function(scrs, centres, conf, col, lty, lwd, ...) {
    mat <- cov.wt(scrs, center = centres)
    if (mat$n.obs == 1)
      mat$cov[] <- 0
    xy <- if (mat$n.obs > 1) {
      veganCovEllipse(mat$cov, mat$center, conf)
    } else {
      scrs
    }
    ordiArgAbsorber(xy, FUN = lines, col = col, lty = lty, lwd = lwd, ...)
  }
  if(missing(main))
    main <- deparse(substitute(x))
  if(missing(sub))
    sub <- paste("method = \"", attr(x, "method"), "\"", sep = "")
  if(missing(xlab))
    xlab <- paste("PCoA", axes[1])
  if(missing(ylab))
    ylab <- paste("PCoA", axes[2])
  t <- if (missing(conf)) {
    1
  } else {
    sqrt(qchisq(conf, df = 2))
  }
  g <- scores(x, choices = axes)
  ng <- length(levels(x$group))
  lev <- levels(x$group)
  ## sort out colour vector if none supplied
  if (is.null(col)) {
    col <- palette()
  }
  col <- rep_len(col, ng)        # make sure there are enough colors
  seg.col <- rep_len(seg.col, ng)     # ditto for segments
  plot(g$sites, asp = 1, type = "n", axes = FALSE, ann = FALSE, ...)
  ## if more than 1 group level
  if(is.matrix(g$centroids)) {
    for(i in seq_along(lev)) {
      curlev <- lev[i]
      take <- x$group == curlev
      j <- which(lev == curlev)
      if (segments) {
        segments(g$centroids[j, 1L], g$centroids[j, 2L],
                 g$sites[take, 1L],
                 g$sites[take, 2L], col = seg.col[i], lty = seg.lty,
                 lwd = seg.lwd)
      }
      if(hull) {
        ch <- chull(g$sites[take, ])
        ch <- c(ch, ch[1])
        lines(x$vectors[take, axes][ch, ], col = col[i], lty = lty,
              lwd = lwd, ...)
      }
      if (ellipse) {
        Ellipse(g$sites[take, , drop = FALSE],
                centres = g$centroids[j, ],
                conf = t,
                col = col[i], lty = lty, lwd = lwd, ...)
      }
      points(g$centroids[j, , drop = FALSE], pch = 16, cex = 1,
             col = col[i], ...)
    }
  } else {
    ## single group
    if (segments) {
      segments(g$centroids[, 1L], g$centroids[, 2L],
               g$sites[, 1L], g$sites[, 2L], col = seg.col,
               lty = seg.lty, ...)
    }
    if(hull) {
      ch <- chull(g$sites)
      ch <- c(ch, ch[1])
      lines(x$vectors[, axes][ch, ], col = col[1L], lty = lty,
            lwd = lwd, ...)
    }
    if (ellipse) {
      Ellipse(g$sites,
              centres = g$centroids,
              conf = t,
              col = col[1L], lty = lty, lwd = lwd,...)
    }
    points(g$centroids[, 1L], g$centroids[, 2L],
           pch = 16, cex = 1, col = col[1L], ...)
  }
  points(g$sites, pch = pch[x$group], cex = cex, col = col[x$group], ...)
  if (label) {
    ordilabel(x, display = "centroids", choices = axes, cex = label.cex)
  }
  localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, ...)
  localAxis(1, ...)
  localAxis(2, ...)
  localBox(...)
  class(g) <- "ordiplot"
  invisible(g)
}

`veganCovEllipse` <-
  function(cov, center = c(0,0), scale = 1, npoints = 100)
  {
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o]))
  }

`ordiArgAbsorber` <- function(..., shrink, origin, scaling, triangular,
                              display, choices, const, truemean, FUN){
  match.fun(FUN)(...) 
}

eqarrowPlot <- function(graph, layout, edge.lty=rep(1, ecount(graph)),
                        edge.arrow.size=rep(1, ecount(graph)),
                        vertex.shape="circle",
                        edge.curved=autocurve.edges(graph), ...) {
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape="none", vertex.label=NA)
  for (e in seq_len(ecount(graph))) {
    graph2 <- delete.edges(graph, E(graph)[(1:ecount(graph))[-e]])
    plot(graph2, edge.lty=edge.lty[e], edge.arrow.size=edge.arrow.size[e],
         edge.curved=edge.curved[e], layout=layout, vertex.shape="none",
         add=TRUE, vertex.label=NA, ...)
  }
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape=vertex.shape, vertex.label=NA, add=TRUE, ...)
  invisible(NULL)
}


mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size/2,size), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}

layout_bernat <- function(g, seed){
  
  l <- layout.circle(g)
  # set.seed(seed)
  y <- sort(l[,2], decreasing = T, index.return=T)$ix
  # y <- sort(l[,2]+runif(length(l[,2]), 0, (max(l[,2])-min(l[,2]))/5), decreasing = T, index.return=T)$ix
  z <- sort(y, decreasing = F, index.return=T)$ix
  ind <- degree(g, mode="in")
  oud <- degree(g, mode="out")
  ord <- sort(ind-oud, decreasing = T, index.return=T)$ix
  
  xbest <-ord[z]
  E <- 1000
  
  for(k in 1:1000){
    a <- ind-oud
    # a[abs(a)<=1] <- 0
    ord <- sort(a+runif(length(a), 0, 0.5), decreasing = T, index.return=T)$ix
    x <- ord[z]
    l <- layout.circle(g, order=x)
    Enew <- cost_intersect(g,l)
    if(Enew<E){
      xbest <- x
      E <- Enew
    }
  }
  l <- layout.circle(g, order=xbest)
  
  return(l)
}

cost_intersect <- function(g, l){
  rownames(l) <- names(degree(g, mode="out"))
  
  lines <- as_edgelist(g)
  lines <- lapply(1:nrow(lines), function(x){
    x1 <- l[lines[x,1], 1]
    y1 <- l[lines[x,1], 2]
    x2 <- l[lines[x,2], 1]
    y2 <- l[lines[x,2], 2]
    
    s <- 0.01
    
    if(x1<x2){
      x1 <- x1+(x2-x1)*s
      x2 <- x2-(x2-x1)*s
    }else{
      x2 <- x2+(x1-x2)*s
      x1 <- x1-(x1-x2)*s
    }
    if(y1<y2){
      y1 <- y1+(y2-y1)*s
      y2 <- y2-(y2-y1)*s
    }else{
      y2 <- y2+(y1-y2)*s
      y1 <- y1-(y1-y2)*s
    }
    
    return(rbind(c(x1,y1),c(x2,y2)))
  })
  
  com <- combn(1:length(lines),2)
  suma <- 0
  for(i in 1:ncol(com)){
    if(segm_intersect(lines[[com[1,i]]],lines[[com[2,i]]])){
      suma <- suma+1      
    }
  }
  return(suma)
}

layout_minover <- function(g, seed){
  E <- 10
  xbest <- 1:length(V(g))
  for(k in 1:100){
    x <- sample(length(V(g)))  
    l <- layout.circle(g, order=x)
    Enew <- cost_intersect(g,l)
    if(Enew<E){
      xbest <- x
      E <- Enew
    }
  }
  l <- layout.circle(g, order=xbest)
  
  return(l)
}
