%% qlanDerivatives.m
% Uses the symbolic toolbox to calculate QLAN model analytic optima under
% the following conditions:
%   1. All channel unitaries are unity
%   2. Links are isoalted, i.e. cross-talk is disregarded
%
% This script calculates:
%   1. Derivatives of lambda, accidentals A, coincidences C, fidelity F,
%   and entangled bit rate RE
%   2. Analytic location and magnitude of the fidelity maximum
%   3. Numeric location and magnitude of the entangled bit rate maximum
%
% Plots generated:
%   1. Plot of RE vs diff(RE)
%   2. Fidelity and normalized EBR using detector values from
%   arxiv_users.json

% clear
syms u lambda A C D etaA etaB dA dB tau RA RB

[~,eff_N,d_N] = loadUsersFromFile('arxiv_users.json');


% Concrete values for plotting
sub_vars = {'etaA', 'etaB', 'dA','dB', 'tau','D'};
sub_vals = [eff_N(1) eff_N(2) d_N(1) d_N(2) 5e-9 4];
domain = [0 3e8];
mysubs = @(x)(subs(x,sub_vars,sub_vals));

% value assumptions
assume(u > 0)
assume(etaA > 0 & etaA <= 1)
assume(etaB>0 & etaB<=1)
assume(lambda >=0 & lambda<=1)
assume(tau > 0)
assume(D > 0)
assume(C > 0)
assume(A > 0)

% Functions
RA = etaA*u/2+dA;
RB = etaB*u/2+dB;
A = 8*tau*(RA)*(RB);
C = etaA*etaB*u;
lambda = 1/(1+A/C);

F = ((D - 1)*lambda + 1)/D;
RE = (A + C)*log2(sqrt(D)*F);
% RE = QLANModel.ebrUniformState();

% Derivatives
dlambda = diff(lambda,u);
dA = diff(A,u);
dC = diff(C,u);
dF = diff(F,u);
dRE = diff(RE,u);

% Get max ebit rate for later normalization
u_maxRE = vpasolve(subs(dRE,sub_vars,sub_vals) == 0, 3e8);
RE_max = eval(subs(RE,[sub_vars 'u'], [sub_vals u_maxRE]));

% Get analytic form of max Fidelity
u_maxF = solve(dF == 0, 'ReturnConditions', true);
F_max = subs(F,'u',u_maxF.u);

%%
% Getting intersection of fid and normalized ebr
RE_norm = mysubs(RE/RE_max);

% Compromise point, i.e. intersection of fidelity and normalized ebr
u_compromise = eval(vpasolve(mysubs(F) == RE_norm,u));

%%
close all
%% Plot Fidelity and EBR
figure()
% yyaxis left
fplot(subs(F,sub_vars, sub_vals),domain);
title("Fidelity & EBR")
xlabel('Link Flux (s^{-1})')
ylim([0 1])
ylabel('Fidelity')
hold
% yyaxis right
ylabel('Fidelity/Normalized EBR')
fplot(subs(RE/RE_max,sub_vars, sub_vals),domain);
temp = ylim;
if max(temp) ~= 0
    temp(1) = 0;
end
ylim(temp)
% set(gca,'XScale', 'log')
legend('Fidelity', 'EBR')

% Showing function parameters on plot
C = [sub_vars;num2cell(sub_vals)];
specifiers = {'%d' '%d' '%d' '%d' '%d' '%.3e'}; %corresponds with sub_vars
txt = "";
for i = 1:size(C,2)
    frmt = ['\n%s=' specifiers{i}];
    txt = strcat(txt, sprintf(frmt, C{:,i}));
end
text(max(domain)*0.8, max(ylim)*0.7,txt)
hold off
%% Plot EBR and dEBR
figure()
xlabel('Link Flux')
yyaxis left
fplot(subs(RE,sub_vars, sub_vals),domain);
temp = ylim;
if max(temp) ~= 0
    temp(1) = 0;
end
ylim(temp)
ylabel('EBR')
hold 
yyaxis right
ylabel('dEBR')
fplot(subs(dRE,sub_vars, sub_vals),domain);
yline(0)

hold off

