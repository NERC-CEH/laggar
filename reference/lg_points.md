# Get number of points within varying distances

Get number of points within varying distances

## Usage

``` r
lg_points(
  x,
  y,
  dist_seq = NULL,
  bound_poly = NULL,
  mindist = NULL,
  maxdist = NULL,
  incdist = NULL
)
```

## Arguments

- x:

  The measurement points

- y:

  The points to be summed

- dist_seq:

  The sequence of distance values

- bound_poly:

  Optional polygon to be kept within

- mindist:

  The minimum distance, ignored if dist_seq is specified

- maxdist:

  The maximum distance, ignored if dist_seq is specified

- incdist:

  The distance increment, ignored if dist_seq is specified

## Value

list containing matrices of the number of points (`buffer_values`),
buffer area (`buffer_area`) and average number of points per unit area
(`value_parea`) for every distance increment. All matrices will have
number of rows equal to the number of measurement points unless
`bound_poly = NULL` in which case `buffer_area` and `points_parea` will
be returned as `NA`.
