#!Rscript --vanilla
#
# exPATT
# Copyright (C) 2024  CEA, CNRS, Paris-Saclay University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cat("PACKAGES INSTALLATION....\n\n")

iPackages <- installed.packages()

# CRAN packages
options(repos = c(CRAN = "http://cran.rstudio.com"))

if(!"lawstat" %in% iPackages) install.packages("lawstat",dependencies=TRUE,)
if(!"gplots" %in% iPackages) install.packages("gplots",dependencies=TRUE)
if(!"reshape2" %in% iPackages) install.packages("reshape2",dependencies=TRUE)
if(!"ggplot2" %in% iPackages) install.packages("ggplot2",dependencies=TRUE)
if(!"pracma" %in% iPackages) install.packages("pracma",dependencies=TRUE)
if(!"RJSONIO" %in% iPackages) install.packages("RJSONIO",dependencies=TRUE)
if(!"gdata" %in% iPackages) install.packages("gdata",dependencies=TRUE)

# BIOCONDUCTOR Packages
options("BioC_mirror"="https://bioconductor.org/")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

if(!"Biostrings" %in% iPackages) BiocManager::install("Biostrings")
if(!"BSgenome" %in% iPackages) BiocManager::install("BSgenome")
if(!"GenomicRanges" %in% iPackages) BiocManager::install("GenomicRanges")
if(!"org.Mm.eg.db" %in% iPackages) BiocManager::install("org.Mm.eg.db")
if(!"org.Hs.eg.db" %in% iPackages) BiocManager::install("org.Hs.eg.db")
if(!"moe430a.db" %in% iPackages) BiocManager::install("moe430a.db")
if(!"moe430b.db" %in% iPackages) BiocManager::install("moe430b.db")

# GITHUB Packages
if(!"devtools" %in% iPackages) install.packages("devtools",dependencies=TRUE)
devtools::install_github("https://github.com/jcaude/exNGS")

cat("DONE.\n")
