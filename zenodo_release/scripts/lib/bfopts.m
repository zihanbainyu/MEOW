function o = bfopts()
% deterministic (quadrature) bayes-factor options; reproducible, no monte carlo
o = bf.options; o.nDimsForMC = 4; o.verbose = false;
end
