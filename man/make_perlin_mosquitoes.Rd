% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{make_perlin_mosquitoes}
\alias{make_perlin_mosquitoes}
\title{Make perlin raster of mosquitoes}
\usage{
make_perlin_mosquitoes(env_dimensions, resolution = 1)
}
\arguments{
\item{env_dimensions}{List of xmn, xmx, ymn, and ymx}

\item{resolution}{Optional, cells will have edge lengths of this size}
}
\description{
Creates a raster to specified dimensions and resolution containing a layer
of perlin noise.
}
