function out = finiteUnitGeoSum(N,r)
%FINITUNITGEOSUM returns a vector of N elements that sum to 1. The
%parameter r determines the decay of the vector values
if ~(abs(r) < 1)
    error("abs(r) must be less than 1")
end
a = (1-r)/(1-r^N);
out = a*r.^(0:(N-1));
end

