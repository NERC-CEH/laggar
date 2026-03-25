# Extract raster values or summaries within concentric buffers ("doughnuts")

This function creates concentric buffer zones (doughnuts) around
measurement points and extracts raster values within each zone. It can
return either raw cell values or summary statistics (e.g., mean, sum)
for each doughnut. Optionally, it can plot the doughnuts for visual
verification. This is highly recommended.

## Usage

``` r
values_within_dist_rast(
  x,
  y,
  dist_seq,
  func = NULL,
  plot_doughnuts = TRUE,
  object_n = 1,
  bg_layer = NULL
)
```

## Arguments

- x:

  `sf` object. Measurement points around which buffers will be created.

- y:

  Raster object (`terra` SpatRaster). The raster to summarize.

- dist_seq:

  Numeric vector. Distances (in map units) for buffering each
  measurement point. Each value represents the outer radius of a
  doughnut. Must be unique and positive.

- func:

  Summary function applied to raster values within each doughnut.
  Default is `NULL`, which returns raw cell values. Can be:

  - A character string naming a function supported by
    [`exactextractr::exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
    (e.g., `"mean"`, `"sum"`, `"median"`).

  - A custom function that takes a numeric vector and returns a single
    value.

- plot_doughnuts:

  Logical. If `TRUE`, plots the doughnuts for visual verification.
  Recommended to ensure buffering behaves as expected.

- object_n:

  Integer. Row number of the measurement point to plot when
  `plot_doughnuts = TRUE`. A single integer or vector of integers.

- bg_layer:

  Background layer to plot along with the objects. Currently only
  accepts `terra::spatRaster`.

## Value

A list with two elements:

- `value`:

  Matrix or vector of extracted values or summary statistics. Rows
  correspond to measurement points; columns correspond to doughnuts.

- `area`:

  Matrix or vector of areas (in square meters) for each doughnut,
  accounting for partial cell coverage.

## Details

The function uses
[`exactextractr::exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
to handle partial cell coverage at polygon edges. When `func = NULL`,
raw cell values are returned along with their corresponding areas. When
`func` is provided, summary statistics are computed for each doughnut.

## See also

[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html),
[`st_buffer`](https://r-spatial.github.io/sf/reference/geos_unary.html),
[`cellSize`](https://rspatial.github.io/terra/reference/cellSize.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: Summarize mean raster values within 100m and 500m doughnuts
result <- values_within_dist_rast(
  x = points_sf,
  y = raster_layer,
  dist_seq = c(100, 500),
  func = "mean",
  plot_doughnuts = TRUE
)

# Access results
result$value   # matrix of mean values per doughnut
result$area    # matrix of areas per doughnut
} # }
```
