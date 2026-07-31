#############################################################################
## Path analysis function
#############################################################################
# Step-by-step path definition along two curves (here, sacral and pubo-ischial 
# curves), starting at (A1,B1) and ending at (An,Bn).
#
# At each step, from the current pair (Ai, Bj), the algorithm looks at the
# up-to-three legal forward moves:
#     - advance on A only:      (i+1, j)
#     - advance on B only:      (i, j+1)
#     - advance on both:        (i+1, j+1)
# and computes the cross-curve landmark distance (Euclidean distance between
# the two landmarks) that each candidate move would produce. It picks the
# move whose resulting pair has the SHORTEST inter-landmark distance -- i.e.
# the move that brings the head of the fetus to engage orthogonally into
# the narrowest space. Assuming the head is presenting the smallest area, this
# is the path it would follow.
#
# Output: a matrix `path` with one row per step and two columns
#         (A_index, B_index), plus the per-step distance and total path
#         length.
#############################################################################

## ------------------------------------------------------------------------
## 1. Distance calculator
## ------------------------------------------------------------------------

euclidean_dist <- function(p, q) sqrt(sum((p - q)^2))

## ------------------------------------------------------------------------
## 2. Alignment
## ------------------------------------------------------------------------
#
# A, B      : (n x 3) matrices of landmark coordinates, in curve order.
# tie_tol   : if two or more candidate moves are within this tolerance of
#             the minimum distance, prefer the diagonal ("advance both")
#             move -- this keeps the path balanced rather than always
#             favouring one curve on exact ties. Set to 0 to disable
#             tie-breaking (first minimum found wins).
#
align_curves <- function(A, B, tie_tol = 1e-9) {
  stopifnot(nrow(A) == nrow(B))
  n <- nrow(A)
  
  i <- 1L; j <- 1L
  path <- matrix(c(i, j), nrow = 1, ncol = 2)
  dists <- euclidean_dist(A[i, ], B[j, ])
  moves_taken <- character(0)
  
  while (!(i == n && j == n)) {
    cand_names <- character(0)
    cand_idx <- list()
    
    if (i < n) { cand_names <- c(cand_names, "A"); cand_idx[["A"]] <- c(i + 1L, j) }
    if (j < n) { cand_names <- c(cand_names, "B"); cand_idx[["B"]] <- c(i, j + 1L) }
    if (i < n && j < n) { cand_names <- c(cand_names, "AB"); cand_idx[["AB"]] <- c(i + 1L, j + 1L) }
    
    cand_d <- sapply(cand_names, function(nm) {
      idx <- cand_idx[[nm]]
      euclidean_dist(A[idx[1], ], B[idx[2], ])
    })
    
    min_d <- min(cand_d)
    tied <- cand_names[abs(cand_d - min_d) <= tie_tol]
    
    if (length(tied) > 1 && "AB" %in% tied) {
      chosen <- "AB"
    } else {
      chosen <- cand_names[which.min(cand_d)]
    }
    
    nxt <- cand_idx[[chosen]]
    i <- nxt[1]; j <- nxt[2]
    path <- rbind(path, c(i, j))
    dists <- c(dists, cand_d[[chosen]])
    moves_taken <- c(moves_taken, chosen)
  }
  
  colnames(path) <- c("A_index", "B_index")
  
  list(
    path = path,
    distances = dists,          # distance at each visited pair, incl. the start
    moves = moves_taken,        # which move ("A","B","AB") was taken at each step
    total_distance = sum(dists)
  )
}


## ------------------------------------------------------------------------
## 3. PCA flattening (for plotting only -- the alignment itself is computed
##    in full 3D by align_curves(), this is purely a visualisation
##    aid)
## ------------------------------------------------------------------------
#
# Projects the pooled A+B landmarks onto the plane of their first two
# principal components. This gives the best 2D "flattening" of the true 3D
# geometry for visualisation purposes
#

flatten_to_2d <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  if (ncol(A) == 2) {
    # already 2D -- nothing to do
    return(list(A2d = A, B2d = B, var_explained = c(1, 1)))
  }
  
  pooled <- rbind(A, B)
  ctr <- colMeans(pooled)
  pooled_c <- sweep(pooled, 2, ctr, "-")
  
  pc <- prcomp(pooled_c, center = FALSE, scale. = FALSE)
  proj <- pooled_c %*% pc$rotation[, 1:2]
  
  n <- nrow(A)
  var_explained <- (pc$sdev[1:2]^2) / sum(pc$sdev^2)
  
  list(
    A2d = proj[1:n, , drop = FALSE],
    B2d = proj[(n + 1):(2 * n), , drop = FALSE],
    var_explained = var_explained
  )
}

## ------------------------------------------------------------------------
## 3. Plotting helper (optional, base graphics)
## ------------------------------------------------------------------------
# A, B here can be the original 3D landmark matrices -- they are PCA-
# flattened to 2D internally before plotting. The alignment result (indices
# and distances) is unaffected; only the visualisation is 2D.

plot_alignment <- function(result, A, B) {
  path <- result$path
  n_steps <- nrow(path)
  
  flat <- flatten_to_2d(A, B)
  A2d <- flat$A2d; B2d <- flat$B2d
  
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))
  
  var_pct <- round(100 * flat$var_explained, 1)
  xlab_txt <- if (ncol(A) > 2) sprintf("PC1 (%.1f%% var)", var_pct[1]) else "dim 1"
  ylab_txt <- if (ncol(A) > 2) sprintf("PC2 (%.1f%% var)", var_pct[2]) else "dim 2"
  
  plot(A2d[, 1], A2d[, 2], type = "l", col = "steelblue", lwd = 2,
       xlim = range(c(A2d[, 1], B2d[, 1])), ylim = range(c(A2d[, 2], B2d[, 2])),
       asp = 1, xlab = xlab_txt, ylab = ylab_txt, main = "Alignment pairs (PCA-flattened to 2D)")
  lines(B2d[, 1], B2d[, 2], col = "firebrick", lwd = 2)
  points(A2d[, 1], A2d[, 2], pch = 16, cex = 0.5, col = "steelblue")
  points(B2d[, 1], B2d[, 2], pch = 16, cex = 0.5, col = "firebrick")
  for (r in seq_len(n_steps)) {
    ai <- path[r, "A_index"]; bj <- path[r, "B_index"]
    segments(A2d[ai, 1], A2d[ai, 2], B2d[bj, 1], B2d[bj, 2], col = "grey60", lty = 2)
  }
  legend("topleft", legend = c("Sacral curve", "Pubo-ischial curve", "Pairing"),
         col = c("steelblue", "firebrick", "grey60"), lty = 1, bty = "n")
  
  plot(seq_len(n_steps), result$distances, type = "b", pch = 16,
       xlab = "Step", ylab = "Cross-curve distance (3D)",
       main = "Inter-landmark distance per step")
}

### Inputting my data

path.species<-list(NULL)
  
for (sp in 1:nlevels(List$Species)){
  dorsal<-resampled.sacral.lm[[sp]] # landmarks on the sacrum for species sp
  ventral<-ventral.lm[[sp]] # landmarks on the ventral end of the canal ellipse for species sp
  path.species[[sp]]<-align_curves(dorsal,ventral)
  png(paste0("Planes_orientations/Steps_",levels(List$Species)[sp],".png"), width = 1100, height = 550,units = "px")
  plot_alignment(path.species[[sp]], dorsal, ventral)
  dev.off()
}
names(path.species)<-levels(List$Species)

## to extract the landmark pairs for the planes, use path.species[[i]]$path, 
# where i is the number of the species needed as in the order of levels(List$Species)
