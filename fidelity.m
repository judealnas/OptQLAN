function f = fidelity(rho,sigma)
%FIDELITY Calculates the fidelity between the density matrices rho and
%   sigma
[sqrtm_rho, res] = sqrtm(rho);
prod_rho_sig = sqrtm_rho * sigma * sqrtm_rho;
[sqrtm_prod, res] = sqrtm(prod_rho_sig);
f = trace(prod_rho_sig)^2;
end

