# Build doughnut geometries from a list of sf data frames

This function takes a list of `sf` data frames representing
progressively buffered geometries and returns a new `sf` object where
each row contains a doughnut-shaped geometry (i.e., the difference
between successive buffers). This

## Usage

``` r
doughnut_builder(poly_list)
```

## Arguments

- poly_list:

  A list of `sf` data frames. Each data frame must have the same number
  of rows, and represent a different buffer level of the same set of
  geometries.

## Value

An `sf` data frame containing doughnut geometries for each input
feature. The first geometry in each row is preserved, and the rest are
doughnut-shaped differences.
