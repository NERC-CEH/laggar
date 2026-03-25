# Plot all the buffered polygons and a single object

Plot all the buffered polygons and a single object

## Usage

``` r
doughnut_checker(
  doughnut_out,
  nbuffers,
  dist_seq,
  nobjects,
  object_n = 1,
  bg_layer = NULL
)
```

## Arguments

- doughnut_out:

  Takes the output from `doughnut_builder`. An `sf` data frame.

- nbuffers:

  The number of buffers implemented.

- dist_seq:

  Numeric vector. Distances (in map units) for buffering each
  measurement point. Each value represents the outer radius of a
  doughnut. Must be unique and positive.

- nobjects:

  The number of measurement points that were originally buffered.

- object_n:

  The measurement point that you want to plot separately. A single or
  vector of integers.

- bg_layer:

  Background layer to plot along with the objects. Currently only
  accepts `terra::spatRaster`.

## Value

Plots showing all buffered objects and the different buffers for the
object chosen in `object_n`. Returns list containing ggplot objects
