# Summarize raster within incremental buffer zones around points

This function calculates raster values within concentric buffer zones
("doughnuts") around a set of measurement points. It optionally applies
a summary function to raster values within each buffer and returns raw
values or aggregated summaries.

## Usage

``` r
lg_rast(
  x,
  y,
  func = NULL,
  dist_seq = NULL,
  mindist = NULL,
  maxdist = NULL,
  incdist = NULL,
  plot_doughnuts = TRUE,
  object_n = 1,
  bg_layer = NULL
)
```

## Arguments

- x:

  An object representing measurement points. Can be an `sf` object or
  convertible to `sf` via internal helper functions.

- y:

  A `SpatRaster` object from the `terra` package representing the raster
  to summarize.

- func:

  Optional. A function or character string specifying how to summarize
  raster values within each buffer. If `NULL`, raw cell values are
  returned. Accepts functions supported by
  [`exactextractr::exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
  or a custom function that takes a numeric vector and returns a single
  value.

- dist_seq:

  Optional numeric vector of buffer distances. If provided, overrides
  `mindist`, `maxdist`, and `incdist`.

- mindist:

  The minimum distance, ignored if dist_seq is specified

- maxdist:

  The maximum distance, ignored if dist_seq is specified

- incdist:

  The distance increment, ignored if dist_seq is specified

- plot_doughnuts:

  Logical. If `TRUE`, plots the doughnut geometries for visual
  verification. Highly recommended for error checking.

- object_n:

  The measurement point that you want to plot separately. A single or
  vector of integers.

- bg_layer:

  Background layer to plot along with the objects. Currently only
  accepts `terra::spatRaster`.

## Value

A list with three elements:

- buffer_values:

  Matrix or vector of raster values (raw or summarized) for each buffer.

- buffer_area:

  Matrix or vector of buffer areas corresponding to each value.

- value_parea:

  Numeric vector of values normalized by buffer area.

## Details

Internally, this function:

1.  Converts `x` to an `sf` object.

2.  Validates that `y` is a `SpatRaster`.

3.  Generates buffer distances using `.index_check()`.

4.  Builds doughnut geometries and extracts raster values using
    [`values_within_dist_rast()`](https://nerc-ceh.github.io/laggar/reference/values_within_dist_rast.md).

5.  Computes value per unit area for each buffer.

## See also

[`values_within_dist_rast`](https://nerc-ceh.github.io/laggar/reference/values_within_dist_rast.md),
[`doughnut_builder`](https://nerc-ceh.github.io/laggar/reference/doughnut_builder.md),
[`doughnut_checker`](https://nerc-ceh.github.io/laggar/reference/doughnut_checker.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: Summarize mean raster values within 3 buffers around points
lg_rast(x = points_sf, y = raster_layer, func = "mean",
        mindist = 100, maxdist = 300, incdist = 100)
} # }
```
