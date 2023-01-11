
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fase

<!-- badges: start -->
<!-- badges: end -->

The goal of fase is to fit a “Functional Adjacency Spectral Embedding”
to an indexed collection of network snapshots. The behavior of each node
is summarized by a latent process, which is fit in a prespecifed spline
basis. The fase package has options to fit either B-spline (with equally
spaced knots) or smoothing spline latent processes.

## Installation

You can install the development version of fase from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("peterwmacd/fase")
```

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(fase) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r} -->
<!-- mean(1:10) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r, echo = FALSE} -->
<!-- plot(1:10) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
