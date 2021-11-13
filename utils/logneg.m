function en = logneg(rho)
%LOGNEG(rho) calculates the log negativity of density matrix rho
en = neg(rho);
en = log2(2*en+1);
end

