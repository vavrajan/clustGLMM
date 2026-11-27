# A Simulated Dataset of Longitudinal Mixed-Type Data

An artificially created dataset by
[`generate_longitudinal_mixed_type_data`](https://vavrajan.github.io/clustGLMM/reference/generate_longitudinal_mixed_type_data.md)
to demonstrate the use of
[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md).
For details about the parameter setting see the example at the bottom of
[`generate_longitudinal_mixed_type_data`](https://vavrajan.github.io/clustGLMM/reference/generate_longitudinal_mixed_type_data.md).

## Usage

``` r
data("longitudinal_mixed_type_data")
```

## Format

A data frame with 400 rows and 32 columns.

- `i`:

  an id variable

- `x`:

  continuous regressor

- `j`:

  order of observation for given id

- `bs1`,`bs2`,`bs3`,`bs4`,`bs5`:

  B-spline basis for `x`

- `f`:

  factor regressor

- `g`:

  a true group label from 1, 2, 3

- `b_num`,`b_poi`,`b_bin`,`b_ord`,`b_cat`:

  random intercept for the outcome of numerical type, analogously for
  other types

- `pred_num`,`pred_poi`,`pred_bin`,`pred_ord`,`pred_cat1`,`pred_cat2`,`pred_cat3`,`pred_cat4`:

  predictor for the outcome of numerical type, analogously for other
  types, categorical has multiple predictors for different categories

- `c_ord1`,`c_ord2`,`c_ord3`,`c_ord4`:

  ordered intercept for level 1, analogously for other levels

- `ynum`,`ypoi`,`ybin`,`yord`,`ycat`:

  an outcome of numerical type, analogously for other types

## Examples

``` r
data(longitudinal_mixed_type_data)
summary(longitudinal_mixed_type_data)
#>        i                x                   j             bs1          
#>  Min.   :  1.00   Min.   :0.0006543   Min.   :1.00   Min.   :0.000000  
#>  1st Qu.: 25.75   1st Qu.:0.2549268   1st Qu.:1.75   1st Qu.:0.000000  
#>  Median : 50.50   Median :0.4745111   Median :2.50   Median :0.004766  
#>  Mean   : 50.50   Mean   :0.4888765   Mean   :2.50   Mean   :0.168857  
#>  3rd Qu.: 75.25   3rd Qu.:0.7492697   3rd Qu.:3.25   3rd Qu.:0.334331  
#>  Max.   :100.00   Max.   :0.9984985   Max.   :4.00   Max.   :0.666642  
#>       bs2                 bs3                 bs4              bs5           
#>  Min.   :0.000e+00   Min.   :0.0000000   Min.   :0.0000   Min.   :0.000e+00  
#>  1st Qu.:2.570e-06   1st Qu.:0.0001401   1st Qu.:0.0000   1st Qu.:0.000e+00  
#>  Median :1.761e-01   Median :0.1099994   Median :0.0000   Median :0.000e+00  
#>  Mean   :2.783e-01   Mean   :0.2323496   Mean   :0.1565   Mean   :8.122e-02  
#>  3rd Qu.:5.574e-01   3rd Qu.:0.4504682   3rd Qu.:0.2946   3rd Qu.:5.300e-07  
#>  Max.   :7.500e-01   Max.   :0.7499922   Max.   :0.6666   Max.   :9.880e-01  
#>  f             g            b_num              b_poi              b_bin        
#>  0:180   Min.   :1.00   Min.   :-1.46314   Min.   :-1.20963   Min.   :-1.2853  
#>  1:220   1st Qu.:1.00   1st Qu.:-0.32161   1st Qu.:-0.15604   1st Qu.:-0.3925  
#>          Median :1.00   Median :-0.01787   Median : 0.09244   Median :-0.1208  
#>          Mean   :1.49   Mean   :-0.01662   Mean   : 0.08911   Mean   :-0.0169  
#>          3rd Qu.:2.00   3rd Qu.: 0.28918   3rd Qu.: 0.42858   3rd Qu.: 0.3041  
#>          Max.   :2.00   Max.   : 0.98786   Max.   : 1.43231   Max.   : 1.7092  
#>      b_ord               b_cat             pred_num          pred_poi       
#>  Min.   :-1.479798   Min.   :-1.06771   Min.   :-2.5579   Min.   :-2.05978  
#>  1st Qu.:-0.290580   1st Qu.:-0.27209   1st Qu.:-1.3022   1st Qu.:-0.45236  
#>  Median : 0.041163   Median : 0.03408   Median :-0.2868   Median : 0.11042  
#>  Mean   : 0.006223   Mean   : 0.02762   Mean   :-0.1670   Mean   : 0.09589  
#>  3rd Qu.: 0.372597   3rd Qu.: 0.35128   3rd Qu.: 1.0693   3rd Qu.: 0.67256  
#>  Max.   : 1.135763   Max.   : 1.22465   Max.   : 2.2541   Max.   : 2.18586  
#>     pred_bin          pred_ord            c_ord1          c_ord2     
#>  Min.   :-1.7337   Min.   :-1.45491   Min.   :-2.50   Min.   :-1.50  
#>  1st Qu.:-0.1493   1st Qu.:-0.06266   1st Qu.:-2.50   1st Qu.:-1.50  
#>  Median : 0.8323   Median : 0.32602   Median :-1.50   Median :-0.50  
#>  Mean   : 0.9409   Mean   : 0.29124   Mean   :-1.99   Mean   :-0.99  
#>  3rd Qu.: 1.9853   3rd Qu.: 0.69869   3rd Qu.:-1.50   3rd Qu.:-0.50  
#>  Max.   : 4.2716   Max.   : 1.45371   Max.   :-1.50   Max.   :-0.50  
#>      c_ord3          c_ord4       pred_cat1   pred_cat2       
#>  Min.   :-0.50   Min.   :0.50   Min.   :0   Min.   :-1.95098  
#>  1st Qu.:-0.50   1st Qu.:0.50   1st Qu.:0   1st Qu.:-0.52434  
#>  Median : 0.50   Median :1.50   Median :0   Median : 0.04830  
#>  Mean   : 0.01   Mean   :1.01   Mean   :0   Mean   : 0.03441  
#>  3rd Qu.: 0.50   3rd Qu.:1.50   3rd Qu.:0   3rd Qu.: 0.65444  
#>  Max.   : 0.50   Max.   :1.50   Max.   :0   Max.   : 1.82601  
#>    pred_cat3           pred_cat4             ynum              ypoi      
#>  Min.   :-1.477781   Min.   :-1.06771   Min.   :-3.3721   Min.   :0.000  
#>  1st Qu.:-0.511948   1st Qu.:-0.27209   1st Qu.:-1.3052   1st Qu.:0.000  
#>  Median : 0.004071   Median : 0.03408   Median :-0.3548   Median :1.000  
#>  Mean   : 0.020832   Mean   : 0.02762   Mean   :-0.1667   Mean   :1.505  
#>  3rd Qu.: 0.519897   3rd Qu.: 0.35128   3rd Qu.: 1.0078   3rd Qu.:2.000  
#>  Max.   : 1.998756   Max.   : 1.22465   Max.   : 4.7516   Max.   :9.000  
#>       ybin        yord    ycat   
#>  Min.   :0.0000   0: 57   0: 77  
#>  1st Qu.:0.0000   1: 55   1:114  
#>  Median :1.0000   2: 59   2:122  
#>  Mean   :0.6375   3: 84   3: 87  
#>  3rd Qu.:1.0000   4:145          
#>  Max.   :1.0000                  
```
