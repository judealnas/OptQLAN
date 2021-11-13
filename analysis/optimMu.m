%% optimMu
% This script solves for the flux parameter that gives the highest
% entangled bit rate in the simple QLAN noise model.

% JML
% 2021.10.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
close all;      clear all;
% 
% f = @(x,y) (1+x)./(4+x).*(4+2*x+4*y)./(1+2*x+4*y) - log(2*exp(1)) ...
%     - log(2)*log2((1+x)./(4+x));
f = @(x,y) (2*x+4*y+1).*log2(1/2*(x.^2+(4*y+4).*x+4*y.^2)./(x.^2+(4*y+1).*x+4*y.^2)) ...
    + 1/log(2)*( (2*x+4*y+4).*(x.^2+(4*y+1).*x+4*y.^2)./(x.^2+(4*y+4).*x+4*y.^2) ...
    - (2*x+4*y+1) );

g = @(x,y) x.^2 + (4*y+1).*x + 4*y.^2;
h = @(x,y) x.^2 + (4*y+4).*x + 4*y.^2;
gg = @(x,y) 2*x + 4*y + 1;
hh = @(x,y) 2*x + 4*y + 4;

dR = @(x,y) gg(x,y).*log2(h(x,y)./(2*g(x,y))) ...
    + 1/log(2) * 1./h(x,y).*(hh(x,y).*g(x,y) - gg(x,y).*h(x,y));
R = @(x,y) g(x,y).*log2(h(x,y)./(2*g(x,y)));
F = @(x,y) h(x,y)./(4*g(x,y));


N = 1000;
x = linspace(0,2,N).';
% dR = f(x,0);

yE = linspace(0,0.5,N).';     % y values to solve for optimum x.
xE = zeros(N,1);
Rm = zeros(N,1);

for j = 1:N
    ff = @(x) f(x,yE(j));
    ff([1e-9 10])
    xE(j) = fzero(ff,[1e-9 10]);
    Rm(j) = R(xE(j),yE(j));
end
%% FIT POLYNOMIAL
ind = yE < 1/4;         % Only perform fitting for y<1/4.
yFit = yE(ind);
xFit = xE(ind);
Rfit = Rm(ind);

pX = polyfit(yFit,xFit,1);
pR = polyfit(yFit,Rfit,1);

xCurve = pX(1)*yE + pX(2);
rCurve = pR(1)*yE + pR(2);

%% PLOTS
% Example curves
y = (0:0.05:0.25).';

figure

subplot(1,2,1)
hold on;    box on;
for j = 1:length(y)
    plot(x,F(x,y(j)),'LineWidth',2)
end
set(gca,'FontName','Arial','FontSize',16)
xlabel('Flux x=\tau_c\mu')
ylabel('Fidelity F')



subplot(1,2,2)
hold on;    box on;
for j = 1:length(y)
    plot(x,max(R(x,y(j)),0),'LineWidth',2)
    legendLabel{j} = ['y=' num2str(y(j))];
end
set(gca,'FontName','Arial','FontSize',16)
xlabel('Flux x=\tau_c\mu')
ylabel('Normalized EBR R=\eta^2R_E/\tau_c')
legend(legendLabel)





figure
subplot(1,2,1)
plot(yE,xE,'LineWidth',2)
hold on;
plot(yE,xCurve,'LineWidth',2,'LineStyle','--','Color',[0.9290 0.6940 0.1250])
set(gca,'FontName','Arial','FontSize',16);
ylabel('Optimal Flux x = \tau_c\mu')
xlabel('Norm. Dark Counts y = \tau_c d/\eta')
hold on
plot(yE,2*yE,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
legend('x_F',['x=' num2str(pX(1)) 'y+' num2str(pX(2))],'x_R')


subplot(1,2,2)
yyaxis left
plot(yE,max(Rm,0),'LineWidth',2);
hold on;
plot(yE,rCurve,'LineWidth',2,'LineStyle','--','Color',[0.9290 0.6940 0.1250])
ylabel('Optimal EBR R = \eta^2 R_E/\tau_c')
ylim([0 0.7])
yyaxis right
plot(yE,(2*yE+1)./(8*yE+1),'LineWidth',2);
ylim([0 1])
ylabel('Optimal Fidelity F')
set(gca,'FontName','Arial','FontSize',16);
% ylabel('Optimal EBR \eta^2 R_E/\tau_c')
xlabel('Norm. Dark Counts y = \tau_c d/\eta')
legend('R_{max}',['R=' num2str(pR(1)) 'y+' num2str(pR(2))],'F_{max}')


