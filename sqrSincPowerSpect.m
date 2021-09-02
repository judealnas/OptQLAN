function h = sqrSincPowerSpect(f0, fwhm, P0)
%h = sqrSincPowerSpect(f0, fwhm, P0) returns function handle h to a squared
%sinc function with center at f0, full-width-half-max of fwhm, and total
%power P0

% half-width of unit normalized sinc^2
% obtained using vpasolve((sin(x)./x).^2 == 0.5, x)
hw = 1.3915573782515101503199613646185; 

% x-axis scaling to obtain desired FWHM
% Additional pi term as MATLAB's sinc(x) is defined as sin(pi*x)/(pi*x)
s = 2*hw/(pi*fwhm); 
h = @(x) P0*s*(sinc(s*(x-f0))).^2;
end

