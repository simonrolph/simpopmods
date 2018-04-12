# R package: simpopmods

An R package for generating simulated population models (IPMs and MPMs)
:herb::mushroom::pig2::seedling::cactus:whale2::deciduous_tree::water_buffalo:

## About

This is an R package I have been developing for my PhD research to simulate population models. It uses integral projection model

## Installation

How to install this package from GitHub

```
install.packages("devtools")
library(devtools)
install_github("simonrolph/simpopmods")
```

## Usage

### Define an IPM

Before you generate any simulated populations, you first need to define a template for your population model. This is done through a IPM descriptor object. You can generate a template for the IPM descriptor by using the `init_IPM_desc()` function.

An `simpopmods` IPM descriptor consists of 
 * Discrete states
 * Continous states (which all share the same domain)
 * Parameters
 * Functions for generating parameters that are not sampled directly
 * Demographic functions (equations of parameters)
 * Kernels (consisting of demographic functions)
 * Functions for the upper limit, lower limit and resolution of the kernel of the continious state domain.
 
### Create an IPM

Create a parameter set or use a sampler to sample parameter sets. I have been using an Adaptive Metropolis Algorithym with delayed rejection from the R package `FME`. This requires some informative priors.

Once you have a parameter set you can create an IPM and create the mega_kernels by:

```
IPM <- make_kernels(params,IPM_desc)
IPM_mega <- make_mega_kernel(IPM)
```

This IPM can be discritised to an MPM by specifying a target stabe stage distribution for the MPM like so:

```
qtiles <- c(0,0.2,0.4,0.7,1)
MPM <- make_MPM (params, IPM_desc, qtiles, submesh_res = 200)
MPM_mega <- make_mega_kernel(MPM)
```

lambda doesn't change in the discretisation process and you can use the `calc_dom_eig()` function to calculate the dominant eigenvalue from the mega kernel:

```
lambda.IPM <- calc_dom_eig(IPM_mega)
lambda.MPM <- calc_dom_eig(MPM_mega)
if (lambda.IPM == lambda.MPM){
  print("it works")
}
```

This R package has tests inplemented with R package `testthat` if you don't believe me.


## Related packages

 * [popbio](https://github.com/cstubben/popbio) - popbio is an R package for modeling population growth rates using age- or stage-classified matrix models. The package consists mainly of methods described in Hal Caswell's Matrix Population Models (2001) and Morris and Doak's Quantitative Conservation Biology (2002).
 * [Rage](https://github.com/jonesor/Rage)
 * [Rcompadre](https://github.com/jonesor/Rcompadre)

## IPM book

[Data-driven Modelling of Structured Populations](http://www.springer.com/gb/book/9783319288918) :star::star::star::star::star: 





