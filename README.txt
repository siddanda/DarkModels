----------------

OVERVIEW AND USAGE

This repository provides plots that are contained in or supplemental to https://arxiv.org/abs/2202.00012. 
If you use any of the information provided in this repository, please cite https://arxiv.org/abs/2202.00012.

----------------

CONTENTS

This repository contains three main folders, with each containing plots for each of the sets
of priors used in the Markov-Chain Monte Carlo sampling runs.

The folder "MCMC Plots" contains plots with information about the specific Markov-Chain
Monte Carlo (MCMC) runs done. "Paper Plots" contains the plots in the main paper at 
https://arxiv.org/abs/2202.00012, as well as related plots with different MCMC priors. 
"Rotation Curves" contain plots of the galactic rotation curves, containing both observed data and
fits via the used models.

Within each of these folders, the subfolders of "Abundance Matching Prior", "Benchmark Analysis", and
"Mass-to-Light Ratio Prior" contain the respective plots to each of those three sets of MCMC priors examined.

In "MCMC Plots", each of these three subfolders contain further subfolders of "DC14", "NFW", and "SIDM",
containing the autocorrelation, corner, and step plots of the MCMC runs done for the respective dark
matter model.

----------------

UNCONVERGED GALAXIES

For each of the MCMC runs, there were a small number of galaxies that did not pass our test 
for autocorrelation. These galaxies should be regarded as exceptions and be treated with extra care in
any analysis. These galaxies are tabulated below:

Abundance Matching Prior
--------
DC14 Unconverged Galaxies - (None)
NFW Unconverged Galaxies - NGC6015
SIDM Unconverged Galaxies - NGC0247

Benchmark Analysis
--------
DC14 Unconverged Galaxies - (None)
NFW Unconverged Galaxies - NGC6015
SIDM Unconverged Galaxies - NGC0247

Mass-to-Light Ratio Prior
--------
DC14 Unconverged Galaxies - NGC3917, UGC08286
NFW Unconverged Galaxies - NGC6015
SIDM Unconverged Galaxies - NGC0289, NGC3521, NGC6674, NGC7841, UGC11914

----------------



