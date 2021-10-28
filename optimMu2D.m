%% optimMu2D.m
% Generate figures plotting the maximum fidelity, maximum EBR and the
% fluxes at which these two optima occur for the dimensionless two
% parameter, single link QLAN model
clear
close all

syms x yA yB D

x_label = 'y_A=\tau_{c}d_A/\eta_A';
y_label = 'y_B=\tau_{c}d_B/\eta_B';
flux_label = 'x=\tau_{c}\mu';
bound_string = '(y_A-y_B)^2-2(y_A+y_B)+1=0';

assume(x >= 0)
assume(yA >= 0)
assume(yB >= 0)
assume(in(D,'integer') & D>1)
D = 4;

h = x^2 + (2*yA + 2*yB + D)*x + 4*yA*yB;
g = x^2 + (2*yA + 2*yB + 1)*x + 4*yA*yB;

F = 1/D*h/g;
REsym = g*log2(sqrt(D)*F);

dF = diff(F,x);
dREsym = diff(REsym,x);

RE = matlabFunction(REsym,'Vars',{x,yA,yB});
dRE = matlabFunction(dREsym,'Vars',{x,yA,yB});

alpha = @(d) (d-sqrt(d))/(sqrt(d)-1);
y_nzRE = @(y,d)y-sqrt(2.*alpha(d).*y)+(alpha(d)./2).^2;

x_Fmax = solve(dF == 0,x);

Fmax = subs(F,x,x_Fmax);


%% Calculating REmax, X for REmax for varying yA, yB and height of bounds
Y = linspace(0,1,500);
[YA, YB] = meshgrid(Y,Y);
REmax = zeros(numel(Y));
x_REmax = zeros(numel(Y));
bound_line = y_nzRE(Y,D); % Points beydon which REmax = 0

for i = 1:numel(Y)
    for j = 1:numel(Y)
        ya = YA(i,j);
        yb = YB(i,j);
        
        % Search for root of dRE; If complex root or NaN returned, search
        % again using different start point
        [xremax] = fzero(@(x)dRE(x,ya,yb),1,optimset('Display','off'));
        xstart = 1/2;
        while imag(xremax) || isnan(xremax)
            xremax = fzero(@(x)dRE(x,ya,yb),xstart);
            xstart = xstart / 2;
        end
        
        remax = RE(xremax, ya,yb);
        
        if remax < 0
            xremax = NaN;
            remax = 0;
        end
        
        x_REmax(i,j) = xremax;
        REmax(i,j) = remax;
    end
end

% Find height of bounding line
bound_height = zeros(size(bound_line));
for i=1:numel(bound_line)
    bound_height(i) = fzero(@(x)dRE(x,Y(i),bound_line(i)),0.5);
end
%% Fitting plane to non-zero REmax
idx_nz = find(REmax > 0);
plane_RE_coeff = [YA(idx_nz), YB(idx_nz) ones(size(idx_nz))]\REmax(idx_nz);
plane_RE = plane_RE_coeff(1).*YA + plane_RE_coeff(2).*YB + plane_RE_coeff(3);
plane_RE(plane_RE < 0) = NaN;
plane_RE_str = compose("%.3fy_A+%.3fy_B+%.3f",plane_RE_coeff(:).');
%% Fitting plane to x_REmax
idx_nnan = idx_nz;
plane_xRE_coeff = [YA(idx_nnan) YB(idx_nnan) ones(size(idx_nnan))]\x_REmax(idx_nnan);
plane_xRE = plane_xRE_coeff(1).*YA + plane_xRE_coeff(2).*YB+plane_xRE_coeff(3);
plane_xRE(plane_xRE < 0) = NaN;
plane_xRE_str = compose("%.3fy_A+%.3fy_B+%.3f",plane_xRE_coeff(:).');
%% Fitting circle to boundary curve
Ycol = Y(:);
boundcol = bound_line(:);
circ_coeff = [Ycol.^2 Ycol ones(numel(Y),1)]\[boundcol.^2 boundcol]
%% Plotting REmax surface, plane fit, and yA == yB line truncating when < 0
fig_REmax = figure('Name','Max EBR');
clf(fig_REmax)
ax_REmax = axes(fig_REmax);

% Construct line of REmax surface where yA == yB and REmax ~=0
y_line = REmax(YA==YB & REmax);
x_line = YA(YA==YB & REmax);

% Perform linear fit
line_coeff = polyfit(x_line, y_line,1);

% Plot maximum RE
surf(ax_REmax,Y,Y,REmax,'EdgeAlpha',0)
hold(ax_REmax, 'on')

% Plot linear fit
plot3(ax_REmax,x_line,x_line,y_line,'LineWidth',1,'Color','k')
linfit_str = compose("x_R=%.2fy+%.2f",line_coeff(:).');
text(min(x_line),min(x_line),max(y_line)-0.02,...
    strcat('$',linfit_str,'$'),'Color','k','Interpreter','latex');

% Plot boundary
plot3(ax_REmax,Y,bound_line,zeros(size(bound_line)),'LineWidth',0.5,'Color','r')
text(ax_REmax,median(Y)+.02, median(bound_line)+0.02,.05,...
    strcat('$',bound_string,'$'), 'Color', 'r','Interpreter','latex','FontWeight','bold')

% Plot planar fit of REmax
% surf(ax_REmax,YA,YB,plane, 'FaceAlpha',0.3, 'EdgeAlpha',0,'FaceColor','r')
hold(ax_REmax, 'off')

xlabel(x_label)
ylabel(y_label)
zlabel('$\mathcal{R}_{max}$','interpreter','latex')
title('$\mathcal{R}_{max}$','interpreter','latex')
view(0,90)
%% x_REmax surface, yA == yB line, and positive EBR lines
fig_xRemax = figure('name','X of REmax');
clf(fig_xRemax)
ax_xRemax = axes(fig_xRemax);

% Line on x_REmax surface where yA == yB, corresponding to simplified 1D
% analysis
y_line = x_REmax(YA==YB & ~isnan(x_REmax));
x_line = YA(YA==YB & ~isnan(x_REmax));

% Linear fit of line
line_coeff = polyfit(x_line,y_line, 1);

% Plot x for max RE
surf(Y,Y,x_REmax,'EdgeAlpha',0)
hold(ax_xRemax,'on')

% Plot planar fit
surf(Y,Y,plane_xRE,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','none')
text(ax_xRemax,0.45,0.85,0,...
    strcat('$',plane_xRE_str,'$'),'Color','b','Interpreter','latex');

% Plot linear fit of yA=yB line
plot3(ax_xRemax,x_line,x_line,y_line,'LineWidth',1,'Color','k')
linfit_str = compose('x_R=%.2fy+%.2f',line_coeff(:).');
text(0.3,0.4,0.8,strcat('$',linfit_str,'$'),'Color','k','Interpreter','latex');

% Plot boundary
plot3(ax_xRemax,Y,bound_line, bound_height,'Color','r','LineWidth',1.5);
text(0.2,0.5,0.6,strcat('$',bound_string,'$'),'Color','r','Interpreter','latex')
hold(ax_xRemax,'off')


title('Optimal $x_\mathcal{R}$','Interpreter','latex')
xlabel(x_label)
ylabel(y_label)
zlabel(flux_label)
view(ax_xRemax,135,45)
%% Fmax surface
fig_Fmax = figure('Name','Max Fidelity');
clf(fig_Fmax)
ax_Fmax = axes(fig_Fmax);

fsurf(ax_Fmax, Fmax, [0 1])

zlim([0 1]);
xlabel('y_A=\tau_{c}d_A/\eta_A')
ylabel('y_B=\tau_{c}d_B/\eta_B')
title('$\mathcal{F}_{max}$','Interpreter', 'latex')

%% x_Fmax surface
fig_xFmax = figure('Name','X of Max Fidelity');
clf(fig_xFmax)
ax_xFmax = axes(fig_xFmax);

fcontour(ax_xFmax, x_Fmax, [0 1])

title("$x_F=2\sqrt{y_Ay_B}$",'Interpreter','latex')
xlabel(x_label)
ylabel(y_label)
zlabel(flux_label)

%% Plotting Xmaxes together
fig_xmax = figure('Name','X Max Fid/EBR');
clf(fig_xmax)
ax_xmax = axes(fig_xmax);

% Plot flux for max fidelity
fsurf(ax_xmax,x_Fmax, [0 1 0 1],'DisplayName','$x_\mathcal{F}$')%,'EdgeColor', 'none');
hold(ax_xmax,'on')

% Plot flux for max EBR
surf(Y,Y,x_REmax,'LineStyle','--','DisplayName','$x_\mathcal{R}$')%,'EdgeColor', 'none');

% Plot bound line
plot3(Y,bound_line,bound_height,'LineWidth',1.5,'Color','r','DisplayName',strcat('$',bound_string,'$'))

xlabel(x_label)
ylabel(y_label)
title('Optimal Flux')
colorbar
legend('Interpreter','latex')
hold(ax_xmax,'off')

%% Plotting max quantities together
fig_maxes = figure();
clf(fig_maxes)
ax_maxes = axes(fig_maxes);

% Calculating boundary of x_REmax

fsurf(ax_maxes, Fmax, [0 1 0 1],'EdgeColor','none');
hold(ax_maxes, 'on')
surf(ax_maxes, Y,Y, REmax,'EdgeColor','none');
% surf(ax_maxes,YA,YB,plane, 'FaceAlpha',0.3, 'EdgeAlpha',0,'FaceColor','r')
plot3(ax_maxes,Y,bound_line,zeros(size(bound_line)),'Color','k','LineWidth',2);
hold(ax_maxes,'off')

%% Plotting RE and F functions gifs
gif_RE = "2D_dimless_RE.gif";
gif_F = "2D_dimless_F.gif";

fig_F = figure('Name', 'Fidelity');
fig_RE = figure('Name','RE');

clf(fig_F);
clf(fig_RE);

ax_F = axes(fig_F);
ax_RE = axes(fig_RE);

F_h = matlabFunction(F,'Vars',{x, yA, yB});

Y1 = linspace(0,1,11);
Y2 = linspace(0,1,11);

domain = [0 2];

for i = 1:numel(Y2)
    yb = Y2(i);
    b_gifoverwrite = i==1;
    for ya = Y1
        fplot(ax_F, @(x)F_h(x, ya, yb), domain, 'DisplayName', sprintf("y_A=%.2f",ya));
        title(ax_F,sprintf("y_B=%.2f",yb))
        
        fplot(ax_RE, @(x)RE(x, ya, yb), domain, 'DisplayName', sprintf("y_A=%.2f",ya));
        title(ax_RE,sprintf("y_B=%.2f",yb))
        
        hold(ax_F,'on')
        hold(ax_RE,'on')
    end
    
    ylim(ax_F, [0.25 1])
    ylim(ax_RE, [0 1])
    
    yticks(ax_F,[0.25 0.3:0.1:1])
    xlabel(ax_F,flux_label)
    xlabel(ax_RE,flux_label)
    
    ylabel(ax_F,'Fidelity')
    ylabel(ax_RE,'Normalized EBR')
    
    legend(ax_F)
    legend(ax_RE)
    
    hold(ax_F,'off')
    hold(ax_RE,'off')
    
    writeToGIF(fig_F,gif_F,b_gifoverwrite, 'DelayTime',0.75)
    writeToGIF(fig_RE,gif_RE,b_gifoverwrite,'DelayTime',0.75)
end

close(fig_F)
close(fig_RE)
%% Plot boundary region
u = linspace(0,1,100);
v = y_nzRE(u,D);
fig_bounds = figure();
clf(fig_bounds);
ax_bounds = axes(fig_bounds);

color = 'r';
area(u,v,'FaceColor',color,'FaceAlpha',0.2,'EdgeColor',color,'LineWidth',1.5,'LineStyle','--');
title("Region where $\mathcal{R}_{max}>0$",'Interpreter','latex')
xlabel("y_A")
ylabel("y_B")
text(0.3,0.5,"$(y_A-y_B)^2-2(y_A+y_B)+1>0$",'Color',color,'Interpreter','latex')


