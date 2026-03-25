# Visualise scalar-on-function regression data

Plots the response, covariate matrix and index matrix for a given SoFR
term

## Usage

``` r
lg_plot(...)

# S3 method for class 'sf'
lg_plot(data, ...)

# S3 method for class 'data.frame'
lg_plot(
  data,
  response = "response",
  covariate = "covariate",
  index = "index",
  widths = c(0.1, 0.5, 0.5),
  ...
)

# S3 method for class 'list'
lg_plot(
  data,
  response = "response",
  covariate = "covariate",
  index = "index",
  widths = c(0.1, 0.5, 0.5),
  ...
)
```

## Arguments

- ...:

  Other arguments passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)

- data:

  Either a data.frame, sf object, or a list

- response:

  string giving the name of the column with the response in it

- covariate:

  string giving the name of the matrix-column with covariate values in
  it

- index:

  string giving the name of the matrix-column with index values in it

- widths:

  the widths of the different plots

## Value

a `patchwork` combination of `ggplot2` objects

## Methods (by class)

- `lg_plot(sf)`: sf interface

- `lg_plot(data.frame)`: data.frame interface

- `lg_plot(list)`: list interface
