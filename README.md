
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FASE

<!-- badges: start -->
<!-- badges: end -->

The goal of fase is to fit a “Functional Adjacency Spectral Embedding”
to an indexed collection of network snapshots. The behavior of each node
is summarized by a latent process, which is fit in a prespecifed spline
basis. The fase package has options to fit either B-spline (with equally
spaced knots) or smoothing spline latent processes. For more details,
see [MacDonald et al., (2022+)](https://arxiv.org/abs/2210.07491).

Additional functions are included in this package to do Procrustes
alignments of latent processes stored as 3-dimensional arrays.

## Installation

You can install fase from CRAN with:

``` r
install.packages("fase")
```

## Example

This is a basic example which shows you how to generate functional
network data and fit a functional adjacency spectral embedding.

First, we generate $50$ snapshots of a functional network on $100$
nodes, with Gaussian edges (standard deviation $\sigma=4$) and
two-dimensional sinusoidal latent processes.

``` r
# load fase package
library(fase)

# set a seed
set.seed(1)

# generate functional network data
data <- gaussian_snapshot_ss(n=100,d=2,
                             x_vec=seq(0,1,length.out=50),
                             self_loops=FALSE,
                             sigma_edge=4)
```

Then we fit a two-dimensional functional adjacency spectral embedding to
this data, using a $B$-spline basis design with dimension $9$. This
choice of basis dimension is chosen to minimize the network generalized
cross validation criterion.

As a post-processing step, and due to the intrinsic non-identifiability
of the latent processes, we align our embedding to the original
processes used to fit the data. Note that by using proc_align3 (rather
than proc_align_slicewise3), we apply the same orthogonal transformation
to all the evaluations of the latent process snapshots. The difference
between slicewise and common identifying transformations is discussed in
more detail in [MacDonald et al.,
(2022+)](https://arxiv.org/abs/2210.07491), Appendix C.

``` r
# embed with fase (B-spline design)
fit <- fase(data$A,d=2,self_loops=FALSE,
            spline_design=list(type='bs',q=9,x_vec=data$spline_design$x_vec))

# align fitted and true latent processes
Z_align <- proc_align3(fit$Z,data$Z)
```

Finally, we extract a couple of example estimated (and aligned) latent
processes: the first dimension for node 34, and the second dimension for
node 54. We plot them against the true values over the entire index
space.

<img src="man/figures/README-plotting-1.png" width="100%" />
