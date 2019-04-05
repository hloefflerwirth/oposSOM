This is the latest development version of the Bioconductor R-package 'oposSOM':

http://master.bioconductor.org/packages/devel/bioc/html/oposSOM.html

To install the package from GitHub use the following commands in R:

```
install.packages("devtools")

library(devtools)
install_github("hloefflerwirth/oposSOM")
```


Note: C++ compiler is required. For Windows systems usage of Rtools is recommended:

https://cran.r-project.org/bin/windows/Rtools/


Alternatively, precompiled package binaries for Windows can be obtained from the latest releases (https://github.com/hloefflerwirth/oposSOM/releases):

```
install.packages(
  "https://github.com/hloefflerwirth/oposSOM/releases/download/2.0.3/oposSOM_2.0.3.zip", 
  repos = NULL, type = "binary" )
```