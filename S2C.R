#' Numerically stable sinc function
#'
#' @param x A numeric vector
#'
#' @return `sin(x)/x`
#' @export
#'
#' @examples
sinc <- function(x) {
  f <- sin(x)
  # Use Taylor expansion close to x=0
  # Minimiser of the relative error:
  # sin(x)/x = 1-x^2/6+x^4/120-O(x^6)...
  # Truncated Taylor: 1-x^2/6
  # Approximate relative computational error for sin(x)/x :
  #   |(1 + eps_1)/(1 + eps_2) - 1| <= 2 * eps
  # Set errors equal (approximatively):
  # x^4 / 120 = 2 * eps
  # x = (240 * eps)^{1/4}
  x0 <- (240 * .Machine$double.eps)^(1/4)
  ok <- abs(x) > x0
  if (all(ok)) {
    f <- f / x
  } else {
    f[ok] <- f[ok] / x[ok]
    f[!ok] <- 1 - x[!ok]^2 / 6
  }
  f
}

#' Title
#'
#' @param omega
#' @param S_fun
#' @param h
#' @param pointwise
#' @param n_folds
#' @param fold If `TRUE`, computes a folded spectrum. If `FALSE`, computes
#' a truncated spectrum.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fold_spectrum <- function(omega, S_fun, h, pointwise = TRUE, n_folds = 10,
                          fold = TRUE,
                          ...) {
  S <- 0
  done <- FALSE
  if (fold) {
    kk <- seq(-n_folds, n_folds, by = 1)
  } else {
    kk <- 0
  }
  k_ <- rep(1, length(h))
  while (!done) {
    omega_ <- omega + kronecker(
      as.matrix(rep(1, NROW(omega))),
      t(2 * pi * kk[k_] / h)
    )
    scaling <- 1
    if (!pointwise) {
      for (d in seq_along(h)) {
        scaling <- scaling * sinc(omega_[, d] * h[d] / 2)^2
      }
    }
    S <- S + S_fun(omega_, ...) * scaling

    idx <- which(k_ < length(kk))
    if (length(idx) == 0) {
      done <- TRUE
    } else {
      idx <- max(idx)
      k_[idx] <- k_[idx] + 1
      if (idx < length(h)) {
        k_[(idx + 1):length(h)] <- 1
      }
    }
  }
  global_scaling <- 1
  for (d in seq_along(h)) {
    global_scaling <- global_scaling *
      (omega[, d] >= -pi / h[d]) *
      (omega[, d] < pi / h[d])
  }
  S <- S * global_scaling
  S
}

#' Title
#'
#' @param x
#' @param dim integer vector of dimension extents, suitable for `array()`.
#' The length of `dim` is the dimension for the intended fft transformation
#'
#' @return
#' @export
#'
#' @examples
fftshift <- function(x, dim = length(x)) {
  stopifnot(prod(dim) == length(x))
  x <- array(as.vector(x), dim = dim)
  if (length(dim) == 1) {
    x <- x[
      c(seq_len(dim[1] / 2) + dim[1] / 2, seq_len(dim[1] / 2)),
      drop = FALSE
    ]
  } else if (length(dim) == 2) {
    x <- x[
      c(seq_len(dim[1] / 2) + dim[1] / 2, seq_len(dim[1] / 2)),
      c(seq_len(dim[2] / 2) + dim[2] / 2, seq_len(dim[2] / 2)),
      drop = FALSE
    ]
  } else {
    stop(paste0("Dimension ", length(dim), " not implemented."))
  }
  x
}
#' Title
#'
#' @param S
#' @param dim
#' @param h
#'
#' @return
#' @export
#'
#' @examples
S2C <- function(S, dim, h) {
  if (length(dim) == 1) {
    fft <- fftwtools::fftw
  } else if (length(dim) == 2) {
    fft <- fftwtools::fftw2d
  } else {
    stop(paste0("Dimension ", length(dim), " is not implemented."))
  }
  C <- fftshift(Re(fft(fftshift(S, dim))), dim)
  C * prod(2 * pi / h / dim)
}

#' Title
#'
#' @param S
#' @param dim
#' @param h
#'
#' @return
#' @export
#'
#' @examples
S2sample <- function(S, dim, h, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(dim) == 1) {
    fft <- fftwtools::fftw
  } else if (length(dim) == 2) {
    fft <- fftwtools::fftw2d
  } else {
    stop(paste0("Dimension ", length(dim), " is not implemented."))
  }
  SD <- sqrt(S * prod(2 * pi / h / dim))
  sample <- Re(fft(fftshift(
    (rnorm(prod(dim)) + 1i*rnorm(prod(dim))) * SD, dim)))
  sample
}



#' Title
#'
#' @param dim
#' @param L
#'
#' @return
#' @export
#'
#' @examples
make_x <- function(dim, L) {
  x <- list()
  for (d in seq_along(dim)) {
    x[[paste0("x", d)]] <-
      (seq_len(dim[d]) - 1 - dim[d] / 2) / dim[d] * L[d]
  }
  x
}

#' Title
#'
#' @param dim
#' @param L
#'
#' @return
#' @export
#'
#' @examples
make_omega <- function(dim, L) {
  h <- L / dim
  w <- list()
  for (d in seq_along(dim)) {
    w[[paste0("w", d)]] <-
      seq(-(dim[d] / 2), dim[d] / 2 - 1, by = 1) / (dim[d] / 2) * pi / h[d]
  }
  w
}





