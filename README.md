# elglm
matlab routines for demonstrating the expected log-likelihood 

**Description:** Computes the maximum "expected log-likelihood" for standard regression 
and linear-non-linear Poisson spiking models. Demonstrates how this estimator can be used 
to approximate the maximum likelihood estimate.   

**Relevant publication :**
[Ramirez,A.D.; Paninski, L., "Fast inference in generalized linear models via expected 
log-likelihoods", *Journal of Computational Neuroscience* (36), 2014]



Dependencies 
==========

) Matlab 

) Statistics and Machine Learning Toolbox

) Mark Schmidt's minFunc toolbox 
https://www.cs.ubc.ca/~schmidtm/Software/minFunc_2007.zip

For the spikeGLMDemo.m

)  GLMspiketools from the Pillow lab 
 https://github.com/pillowlab/GLMspiketools/archive/old_v1.zip

) Download and instal GLM Net 
https://web.stanford.edu/~hastie/glmnet_matlab/glmnet_matlab.zip

Installation
===========
1. Download the elglm zip file or clone the repository: "git clone https://github.com/alxdroR/elglm"
2. Download and install the dependent code (minFunc and optionally GLMspiketools and glmnet). 
   GLMspiketools requires mex file compilation. glmnet might as well depending on your version of 
   Matlab and OS. 
3. Add elglm and dependent code to the Matlab path 

USE 
======
Open and read the demo scripts in /demos/ to see simple examples of code usage. 
LNPdemo - compares the maximum likelihood estimator and maximum expected likelihood estimator 
          with and without an L2 penalty on the likelihood. 

 
