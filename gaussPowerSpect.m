function f = gaussPowerSpect(f0, fwhm, P0)
%GAUSSIANPOWERSPECT returns a function handle to a Gaussian shaped power
%spectrum with the specified center frequency f0, full-width half-max
%fwhm, and total power P0.

sigma = fwhm/(2*sqrt(2*log(2)));
f = @(x) P0/(sigma*sqrt(2*pi))*exp(-(x-f0).^2/(2*sigma^2));
end

