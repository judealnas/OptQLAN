function spect = powerSpectrum(f0, fwhm, P0, shape)
%POWERSPECTRUM returns a handle to square integrable function with the
%specified center, FWHM, total, and shape
%
%   SPECT = powerSpectrum(F0, FWHM, P0, SHAPE) returns a function handle to
%   a square integrable function with center F0, full-width-half-max of
%   FWHM. The integral of this function over the domain is P0. SHAPE can be
%   one of:
%       'gauss' - Gaussian spectrum (default)
%       'uniform' - uniform spectrum
%       'sincsqr' - squared sinc shaped spectrum

p = inputParser;
p.addRequired('f0', @(x)isnumeric(x) && isscalar(x));
p.addRequired('fwhm', @(x)isnumeric(x) && isscalar(x));
p.addRequired('P0',@(x)isnumeric(x) && isscalar(x));
p.addOptional('shape', 'gauss', @(x)(isstring(x) && isscalar(x)) || ischar(x));

if nargin == 3
    p.parse(f0,fwhm,P0)
elseif nargin == 4
    p.parse(f0,fwhm,P0,shape)
end

f0 = p.Results.f0;
fwhm = p.Results.fwhm;
P0 = p.Results.P0;
shape = lower(p.Results.shape);

switch shape
    case 'gauss'
        sigma = fwhm/(2*sqrt(2*log(2)));
        spect = @(x) P0/(sigma*sqrt(2*pi))*exp(-(x-f0).^2/(2*sigma^2));
    case 'uniform'
        fL = f0 - fwhm/2;
        fH = f0 + fwhm/2;
        spect = @(x) (x >= fL & x <= fH).*P0/fwhm;
    case 'sincsqr'
        % half-width of unit normalized sinc^2
        % obtained using vpasolve((sin(x)./x).^2 == 0.5, x)
        hw = 1.3915573782515101503199613646185; 

        % x-axis scaling to obtain desired FWHM
        % Additional pi term as MATLAB's sinc(x) is defined as sin(pi*x)/(pi*x)
        s = 2*hw/(pi*fwhm); 
        spect = @(x) P0*s*(sinc(s*(x-f0))).^2;
end