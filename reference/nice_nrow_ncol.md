# Determine Optimal Matrix Dimensions Given the Number of Cells

Did you ever want to arrange *n* small plots to one big plot of matrix
grid? To do so, the easiest way is to set `par(mfrow = c(nrow, ncol))`.
Once, you can set it manually. But what if you need to do this over and
over for different *n* values? The following function proposes `nrow`
and `ncol` such that `nrow * ncol >= n` and the remainder
`nrow * ncol - n` is the smallest possible.

## Usage

``` r
nice_nrow_ncol(n, increasing = TRUE, maxsteps = 20, maxratio = 3)
```

## Arguments

- n:

  number of cells to be covered.

- increasing:

  a Boolean, if TRUE then the values are returned in increasing
  (non-decreasing) order, if FALSE then `nrow >= ncol`.

- maxsteps:

  how many integers away from `sqrt{n}` to try.

- maxratio:

  what is the maximal ratio between `nrow` and `ncol`?

## Details

Function starts at integer `ceiling{sqrt{n}}` and for each integer `x`
between `ceiling{sqrt{n}}` and `ceiling{sqrt{n}} + maxsteps` finds the
lowest integer `y` such that `x*y >= n`. Among all found pairs of `x`
and `y` that fulfill the condition `x/y <= maxratio` it chooses the one
that has the smallest remainder `x*y - n`. If more pairs achieve the
optimum, the pair with the smallest ratio `x/y` is chosen. Depending on
`increasing` the optimal pair `c(y,x)` or `c(x,y)` is returned.

## Value

Two integers `c(nrow, ncol)`.

## Author

Jan Vávra

## Examples

``` r
nice_nrow_ncol(10)
#> [1] 2 5
nice_nrow_ncol(11)
#> [1] 3 4
nice_nrow_ncol(15)
#> [1] 3 5
nice_nrow_ncol(15, increasing = FALSE)
#> [1] 5 3

# (2, 17) does not fulfill the condition on maxratio
nice_nrow_ncol(34, maxratio = 3)
#> [1] 5 7

# (2, 17) fulfills the condition on maxratio
nice_nrow_ncol(34, maxratio = 10)
#> [1]  2 17

# ceiling(sqrt(34)) = 6 but 6+5=11 < 17 and pair (2, 17) is not even considered
nice_nrow_ncol(34, maxsteps = 5, maxratio = 10)
#> [1] 5 7
nice_nrow_ncol(34, maxsteps = 11, maxratio = 10)
#> [1] 5 7
nice_nrow_ncol(34, maxsteps = 12, maxratio = 10)
#> [1]  2 17
```
