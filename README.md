# sGMRFmix

[![Build Status](https://travis-ci.com/AntixK/sGMRFmix.svg?branch=dev)](https://travis-ci.com/AntixK/sGMRFmix)

Python library for sGMRFmix model for anomaly detection in time-series data.
sGMRFmix is short for sparse mixture of Gaussian Markov Random Fields.
This is essentially a C++ (and python) port of the R package [`sGMRFmix`](https://cran.r-project.org/web/packages/sGMRFmix/index.html) to make it run faster for larger datasets.

## Model Overview
sGMRFmix is a mixture of GMRFs that predict the likelihood of a random variable using the variables in its markov blanket. Lower the log likelihood, higher the anomaly score. The markov blanket is estimated using a Gaussian graphical model with constraint that enforces sparsity in the inverse covariance matrices of the mixture of GMRF model. This can be done in a stright-forward manner using Graphical LASSO model. You can check out the [paper](https://ide-research.net/papers/2016_ICDM_Ide.pdf) for further details and the math.
  
![sGMRFmix Model](https://github.com/AntixK/sGMRFmix/blob/main/assets/model_overview.png)

## Performance Comparison
The follwing plot shows the performance comparison with the only [`sGMRFmix`](https://cran.r-project.org/web/packages/sGMRFmix/index.html) library available publicaly in R language.
 
![sGMRFmix Model](https://github.com/AntixK/sGMRFmix/blob/dev/assets/sgmrf_comparison.png)

## Installation for google colaboratory
### install armadillo-10.2.1
download armadillo-10.2.1 https://osdn.net/frs/g_redir.php?m=jaist&f=arma%2Farmadillo-10.2.1.tar.xz 

```
cd armadillo-10.2.1
make install
```
### install sgmrfmix
```
cd sGMRFmix
python3 setup.py install
```



## Acknowledgements
- T. Ide, A .Khandelwal, J .Kalagnanam, **Sparse Gaussian Markov Random Field Mixtures for Anomaly Detection**, IEEE 16th International Conference on Data Mining (ICDM), 2016, pp 955–960
- https://rdrr.io/cran/sGMRFmix/f/vignettes/sGMRFmix.Rmd
- https://cran.r-project.org/web/packages/sGMRFmix/vignettes/sGMRFmix.html
- https://github.com/cran/sGMRFmix
- https://github.com/JClavel/glassoFast
