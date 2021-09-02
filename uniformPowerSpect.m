function h = uniformPowerSpect(f0, fwhm, P0)
%h = uniformPowerSpect(f0, fwhm, P0) returns a function handle h for a
% uniform power spectrum centered at f0, width of fwhm, and an amplitude of
% P0/fwhm
fL = f0 - fwhm/2;
fH = f0 + fwhm/2;
h = @(x) (x >= fL & x <= fH).*P0/fwhm;
end

