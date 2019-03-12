Class 7: Functions and Packages
================
Karen Li
1/29/2019

Revisiting Functions
--------------------

``` r
source("http://tinyurl.com/rescale-R")
```

Let's try the rescale() function out

``` r
rescale(c(1, 5, 10))
```

    ## [1] 0.0000000 0.4444444 1.0000000

Let's try **rescale2()** with the **stop()** function catch for non-numeric input

``` r
rescale2(c(1:5), "string")
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
x <- c(3, 7, NA, 4, 8, NA)
which(is.na(x))
```

    ## [1] 3 6

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
sum(is.na(x))
```

    ## [1] 2

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum( is.na(x) & is.na(y) )
```

    ## [1] 1

Take my working snippet and make my first function

``` r
both_na <- function(x, y) {
 sum( is.na(x) & is.na(y) )
}

both_na(x, y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
y3 <- c( 1, NA, NA, NA, NA)
y4 <- c( 1, NA, NA, NA, NA, NA)

both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na(x, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
both_na(x, y4)
```

    ## [1] 5
