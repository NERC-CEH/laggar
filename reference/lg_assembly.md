# Assembling

Assembling

## Usage

``` r
lg_assembly(
  response,
  covariate,
  index_seq = NULL,
  minindex = NULL,
  maxindex = NULL,
  incindex = NULL
)
```

## Arguments

- response:

  Vector of response values

- covariate:

  Matrix of covariate values, assembled by lg\_ summary functions

- index_seq:

  Sequence of index values used to index the covariate matrix

- minindex:

  Minimum index, ignored if `index_seq` specified

- maxindex:

  Maximum index, ignored if `index_seq` specified

- incindex:

  index increment, ignored if `index_seq` specified

## Value

list
