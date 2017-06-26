#' PACKAGE CREATION IN RSTUDIO
#'
#' Install RStudio and GIT
#' Create a GITHUB account
#' In Rstudio/Tools/Global Option set git.exe adress
#' Generate Secure Key and copy it on your github account
#' In Rstudio/Tools/Project Option/Git, set your control system to git
#' Install packages "devtools" and "roxygen2"
#'
#' Create Rstudio new project/package
#' Write your function in ./R/xxxxx_file.R
#' use "document()" to compile your package
#'
#' To sync with github :
#' Create a new folder on github named as your package
#' On your local package folder : right click and open git shell
#' Follow instructions on site to link your local folder to github
#'
#' In Rstudio you can now push or pull your commit (version/revision)
#' in the Git window (right). And you can check/load your package in
#' the Build window.
#'
#' To install the package :
#' devtools::install_github("Mystilivia/SDjoygret")

```{r}
require(devtools) ; require(roxygen2)
devtools::document()
```

```{r}
install_github("Mystilivia/SDjoygret")
require(SDjoygret)
```

