
# matrixNormal 0.0.5 2021-Apr-1

  - Minor Change: Changed lazy data to FALSE to remove notes on CRAN
    check.

# matrixNormal 0.0.4 2020-08-26 on CRAN

  - New `vech()` function: performs half-vectorization on a symmetric
    matrix. This is the vector-form of the lower triangular matrix only.
    Unlike other functions on CRAN, vech() inherits any names from the
    matrix.
  - `dmatrixnorm`: Clarified by replacing the name of argument `use.log`
    with `log` for consistency in argument name with mvtnorm and stats
    package.
  - `is.symmetric()`, `is.positive.definite()`,
    `is.positive.semi.definite()`:
      - if A is not symmetric, these functions NOW return FALSE instead
        of stopping the function. Restructured helper `find.eval()`.  
      - if A contains a missing value (`NA`), these functions NOW return
        NA.
  - `rmatrixnorm`: Added the first argument `s` to draw many random
    samples. Only 1 sample is still drawn; the argument currently has no
    effect but acts as a placeholder in future releases.
  - Clarified documentation.
  - Added session information and version details to the vignette. The
    updated versions of the packages listed do not affect the
    conclusions.

# matrixNormal 0.0.2 2019-12-5

  - The documentation is clarified.
  - Submitted to CRAN but failed to be released.

# matrixNormal 0.0.1 2019-07-09

  - Vignette is added, replacing the package R file. Included the use
    and uniqueness of matrixNormal distribution. Also included an
    example that uses the package. Added additional packages to be
    imported in Documentation.
  - Documentation clarified.
  - Minor Changes in Examples: \*\* Changed order of examples in
    rmatnorm() for reproducibility \*\* Made clearer in using the
    dataset package in pmatnorm() example \*\* Removed installation of
    matrixcalc package in is.symmetric.matrix examples.
  - `is.positive.definite()`, `is.positive.semi.definite()` returns NA
    if the matrix contains missing value (bug fix).
  - `pmatnorm()`, `dmatnorm()`, `rmatnorm()` now throws error if the
    parameters of the matrix Normal Distribution `M`, `U`, or `V`
    contain any missing values.
  - `rmatnorm()` now returns a matrix with rownames from U and the
    colnames from V.
  - In `rmatnorm()`, `pre0.9_9994` that was passed from `rmvnorm()`
    function in **mvtnorm library** is removed, because it is not
    needed. This argument was introduced in mvtnorm library to fix a bug
    in version 0.9-9993, but matrixNormal uses a version of at least
    1.0.8. This argument is just not needed, and if pre0.9\_9994 is set
    to TRUE, nothing will happen.

# matrixNormal 0.0.0.9000

  - This is a new submission.
  - Added a `NEWS.md` file to track changes to the package.
  - First Release of the Package
  - Successfully passed windows check.
