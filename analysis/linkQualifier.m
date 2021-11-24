%% linkQualifier.m
% This script shows that the max EBR of a link depends only on the value of
% (yA-yB)^2 -2(yA+yB)+1 when D=4.
%
% This is shown by finding (yA,yB) pairs that satisfy
% (yA-yB)^2-2(yA+yB)+1=C for arbitrary values of C and plotting the maximum
% EBR values
clear; clc ; close all
%% Generating yA,yB values
D = 4;
Cs = linspace(0,1,6);
func_yb = arrayfun(@(c)@(yA)yA+1-sqrt(4.*yA+c),Cs,'UniformOutput',false);

yAs = linspace(0,1,70);
yBs = arrayfun(@(i)func_yb{i}(yAs),1:numel(func_yb),'UniformOutput',false);
yBs = cell2mat(yBs.'); % row index = C val; col index = yA val
%% Calculating EBRmax values and xEBRmax values
ebrmaxes = zeros(size(yBs));
ebrXmaxes = zeros(size(yBs));
for i = 1:size(yBs,1)
    for j = 1:numel(yAs)
        yA = yAs(j);
        yB = yBs(i,j);
        
        if yB < 0
            ebrmaxes(i,j) = NaN;
            ebrXmaxes(i,j) = NaN;
        else
            model = QLANLinkModel([yAs(j), yBs(i,j)],D);
            [ebrmaxes(i,j), ebrXmaxes(i,j)] = model.EBRMaximum;
        end
        
    end
end
%% Calculate mean and std for each value of C
ebrmax_mean = mean(ebrmaxes,2,'omitnan');
ebrmax_std = std(ebrmaxes,0,2,'omitnan');
log10std = log10(ebrmax_std);
mean_std_str = compose("R_{Emax}=%.2f\\pm%.2E",[ebrmax_mean(:) ebrmax_std(:)]);

fig_max = figure(1);
plot3(yAs, yBs, ebrmaxes);
xlabel("$y_A$")
ylabel("$y_B$")
legend(compose('Q=%.2f ; %s',[num2cell(Cs(:)),mean_std_str]),'Location','best','Interpreter','tex');
view(0,90)
title("$\mathcal{R}_{Emax}$",'Interpreter','latex')
