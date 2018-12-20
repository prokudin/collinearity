# collinearity
We analyze HERMES multiplicities in semi-inclusive deep-inelastic scattering and extract non-perturbative transverse momentum dependence of unpolarized transverse momentum dependent distributions. We discuss the importance of data selection in conducting the fits of the multiplicities.  In particular   we implement   the  collinearity criteria introduced in reference~\cite{Boglione:2016bph} that allow us to study the effects of selecting data   that is  predominantly in the current fragmentation region.  We compare our parameters to previous extractions in order to better interpret our results.  We also give an outlook on what impact this criteria can have for on-going and future experiments.


If you want to produce MC results on your own (using different cuts if need be)
do the following:

source setup.bash, or source setup.csh

jam3d -t 1 inputs/(inputname).py
  
(see list of awailable inputs in the directory ./input)
  
This will produce a single fit, so that you know what parameers you expect.
The input file will be automatically rewritten with new parameters.

Proceed with Monte Carlo runs

jam3d -t 3 -p inputs/(inputname).py
  
mc runs will be created in the directory mcdata

cd mcdata

concatenate all runs:

mcproc .

mcproc -c 0.001

cd ..

copy all final runs into the desired directory:

cp ../(final name).mcp ./samples/(yoour desired final name).mcp
cp ../(summary name).mcp ./samples/(yoour desired summary name>).mcp  
  
clean up

rm ./mcdata/*.*

pkill python

  
