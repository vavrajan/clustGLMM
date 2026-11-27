# Slategray Color Palette

Creates a vector of `n` slategray colors.

## Usage

``` r
slategray.colors(n, start = 0.3, end = 0.9)
slategrey.colors(n, start = 0.3, end = 0.9)
```

## Arguments

- n:

  the number of slategray colors to be in the palette.

- start:

  starting slategray level in the palette (should be between 0 and 1
  where zero indicates "slategray4" and one indicates "slategray1").

- end:

  ending gray level in the palette.

## Details

The function `slategray.colors` chooses a series of `n` slategray levels
between `start` and `end` by interpolating 0=`slategray4`,
1/3=`slategray3`, 2/3=`slategray2`, `slategray1` = 1. The returned
palette contains the corresponding gray colors. This palette is used in
[`plot_cat_vs_x_grouped`](https://vavrajan.github.io/clustGLMM/reference/plot_cat_vs_x_grouped.md)
for levels of group 0.

`slategrey.colors` is an alias for `slategray.colors`.

## Value

A vector of `n` slategray colors.

## Author

Jan Vávra

## Examples

``` r
slategray.colors(5, start = 0.3, end = 0.9)
#> [1] "#98B0C6" "#A6C0D8" "#B1CDE7" "#BAD6F2" "#C0DDF9"
slategrey.colors(5, start = 0.1, end = 0.7)
#> [1] "#7A8C9E" "#90A7BC" "#A2BBD3" "#ADC8E2" "#B8D4EF"
```
