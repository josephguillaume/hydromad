<!-- this document has been modified from the rOpenSci and dplyr CONTRIBUTING files --->
<!-- https://github.com/ropensci/dotgithubfiles/blob/master/dotgithub/CONTRIBUTING.md -->
<!-- https://github.com/tidyverse/dplyr/blob/master/.github/CONTRIBUTING.md -->

# Contributing to hydromad

This document outlines how to contribute changes to hydromad.

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before attempting to make a change. 

Please note we have a code of conduct that we request you adhere to in all your interactions with this project (see below).

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue using our `ISSUE_TEMPLATE.md` and provide steps on how to reproduce the bug.

Prior to making a pull request, we highly recommend you refer to `PULL_REQUEST_TEMPLATE.md` for a full checklist of tasks that needs to be completed.

### Pull request process

The pull request process for contributing to hydromad is as follows:

1. Create a branch in your Git fork and label it appropriately.
2. Open a pull request to hydromad's `master` branch.
3. Fill in the pull request template as required. The pull request template will automatically appear once the request has been made.
4. Discuss any amendments your pull request requires with owners of the repository that need to be made before your pull request can be approved.

#### Other notes

*  We recommend that R code is in the tidyverse [style guide](http://style.tidyverse.org), with the exception of naming using snake_case. Naming conventions are being formalised - in the mean time, please use existing code as an example.
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles to your code. C/C++ should be formatted using [clang-format](https://github.com/llvm-mirror/clang/tree/master/tools/clang-format). **Please don't restyle code that has nothing to do with your pull request**.  
*  We recomend code linting using [lintr](https://github.com/jimhester/lintr) or [clang-tidy](https://github.com/llvm-mirror/clang-tools-extra/tree/master/clang-tidy).
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2) for R documentation - this applies to all new contributions. Legacy code is being updated to conform.
*  We use [testthat](https://cran.r-project.org/package=testthat) for R unit testing. Contributions
with test cases included are easier to accept.  
*  For user-facing changes where the pull request has been approved, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/pull request(s).

#### See rOpenSci [development guide](https://devguide.ropensci.org/) and [CRAN policies](https://cran.r-project.org/web/packages/policies.html) for further details.

## Code of Conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms. Note that `hydromad` is not yet part of ROpenSci, such that not all reporting procedures are yet available.
