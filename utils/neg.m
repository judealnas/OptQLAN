function z = neg(rho)
    rhoPT = zeros(size(rho));   % Initialize partial transpose.
    d = sqrt(size(rho,1));
    for r=1:d:d*(d-1)+1         % Cycle through rows.
        for c=1:d:d*(d-1)+1     % Columns.
            rhoPT(r:r+d-1,c:c+d-1) = rho(r:r+d-1,c:c+d-1).';
        end
    end
    lam = real(eig(rhoPT));  % Compute eigenvalues; remove residual imaginary component.
    z = -(lam<0).'*lam;   % Add up absolute values of negative eigenvalues.
end