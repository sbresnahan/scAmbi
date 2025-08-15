## ----echo=F,message=F,warning=F-----------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 7,  # Width of the plot in inches
  fig.height = 5  # Height of the plot in inches
)

## ----eval=F-------------------------------------------------------------------
#  install.packages(c("knitr", "rmarkdown", "testthat"))

## ----eval=F-------------------------------------------------------------------
#  # install.packages("remotes")
#  remotes::install_github("sbresnahan/scAmbi")

## ----eval=F-------------------------------------------------------------------
#  # If you have scAmbi.zip or a source tar.gz:
#  install.packages("scAmbi.zip", repos = NULL, type = "source")
#  # or use devtools:
#  # devtools::install("path/to/scAmbi/")

## ----include-plotA, echo=FALSE, fig.cap="Figure 1", out.width="100%"----------
knitr::include_graphics("plotA.png")

## ----include-plotB, echo=FALSE, fig.cap="Figure 2", out.width="100%"----------
knitr::include_graphics("plotB.png")

## ----include-plot1, echo=FALSE, fig.cap="Figure 3", out.width="100%"----------
knitr::include_graphics("plot1.png")

## ----include-plot2, echo=FALSE, fig.cap="Figure 4", out.width="100%"----------
knitr::include_graphics("plot2.png")

## ----include-plot3, echo=FALSE, fig.cap="Figure 5", out.width="100%"----------
knitr::include_graphics("plot3.png")

