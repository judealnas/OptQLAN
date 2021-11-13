%% optimMu2D.m
% Generate figures plotting the maximum fidelity, maximum EBR and the
% fluxes at which these two optima occur for the dimensionless two
% parameter, single link QLAN model
clear
close all
set(0,'defaultTextInterpreter','late');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

syms x yA yB D

y1_label = '$y_1$';
y2_label = '$y_2$';
flux_label = '$x=\tau_{c}\mu$';
bound_string = '$(y_1-y_2)^2-2(y_1+y_2)+1=0$';

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

alpha = @(d) (d-sqrt(d))/(sqrt(d)-1)/2;
y_nzRE = @(y,d,c)y+alpha(d)-sqrt(4.*alpha(d).*y+c);

x_Fmax = solve(dF == 0,x);

Fmax = subs(F,x,x_Fmax);


%% Calculating REmax, X for REmax for varying yA, yB and height of bounds
Y = linspace(0,1,500);
[YA, YB] = meshgrid(Y,Y);
REmax = zeros(numel(Y));
x_REmax = zeros(numel(Y));
bound_line = y_nzRE(Y,D,0); % Points beydon which REmax = 0

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
plane_RE_str = compose("%.3fy_1+%.3fy_2+%.3f",plane_RE_coeff(:).');
%% Fitting plane to x_REmax
idx_nnan = idx_nz;
plane_xRE_coeff = [YA(idx_nnan) YB(idx_nnan) ones(size(idx_nnan))]\x_REmax(idx_nnan);
plane_xRE = plane_xRE_coeff(1).*YA + plane_xRE_coeff(2).*YB+plane_xRE_coeff(3);
plane_xRE(plane_xRE < 0) = NaN;
plane_xRE_str = compose("%.3fy_1+%.3fy_2+%.3f",plane_xRE_coeff(:).');
%% Fitting circle to boundary curve
Ycol = Y(:);
boundcol = bound_line(:);
circ_coeff = [Ycol.^2 Ycol ones(numel(Y),1)]\[boundcol.^2 boundcol];
%% Plotting REmax surface, and yA == yB line truncating when < 0
fig_REmax = figure('Name','Max EBR');
clf(fig_REmax)
ax_REmax = axes(fig_REmax);

% Construct line of REmax surface where yA == yB and REmax ~=0
y_line = REmax(YA==YB & REmax);
x_line = YA(YA==YB & REmax);

% Perform linear fit
line_coeff = polyfit(x_line, y_line,1);

% Plot maximum RE
surf(ax_REmax,YA,YB,REmax,'EdgeAlpha',0)
hold(ax_REmax, 'on')
colorbar('FontName','Times New Roman')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
% Plot linear fit
% plot3(ax_REmax,x_line,x_line,y_line,'LineWidth',1,'Color','k')
% linfit_str = compose("x_R=%.2fy+%.2f",line_coeff(:).');
% text(min(x_line),min(x_line),max(y_line)-0.02,...
%     strcat('$',linfit_str,'$'),'Color','k','Interpreter','latex');

% Plot boundary
plot3(ax_REmax,Y,bound_line,zeros(size(bound_line)),'LineWidth',1,'Color','r')
[~,ind] = min(abs(Y-0.2),[],'all','linear');
text_x = Y(ind);
text_y = bound_line(ind);
annotation('textarrow',[0.47,0.37],[0.548,0.269],'String',...
    strcat('$',bound_string,'$'), 'Color', 'r','Interpreter','latex'...
    ,'FontWeight','bold'...#
    ,'LineWidth',1 ...
    ,'FontSize',18)

% Plot planar fit of REmax
% surf(ax_REmax,YA,YB,plane, 'FaceAlpha',0.3, 'EdgeAlpha',0,'FaceColor','r')
hold(ax_REmax, 'off')

xlabel(y1_label)
ylabel(y2_label)
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
colorbar('FontName','Times New Roman')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

grid(ax_xRemax,'off');
% Plot planar fit
% surf(Y,Y,plane_xRE,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','none')
% text(ax_xRemax,0.45,0.85,0,...
%     strcat('$',plane_xRE_str,'$'),'Color','b','Interpreter','latex');

% Plot linear fit of yA=yB line
% plot3(ax_xRemax,x_line,x_line,y_line,'LineWidth',1,'Color','k')
% linfit_str = compose('x_R=%.2fy+%.2f',line_coeff(:).');
% text(0.3,0.4,0.8,strcat('$',linfit_str,'$'),'Color','k','Interpreter','latex');

% Plot boundary
plot3(ax_xRemax,Y,bound_line, bound_height,'Color','r','LineWidth',1.5);
% text(0.2,0.5,0.6,strcat('$',bound_string,'$'),'Color','r','Interpreter','latex')
annotation('textarrow',[0.47,0.37],[0.548,0.269],'String',...
    strcat('$',bound_string,'$'), 'Color', 'r','Interpreter','latex'...
    ,'FontWeight','bold'...
    ,'LineWidth',1 ...
    ,'FontSize',18)

hold(ax_xRemax,'off')

title('$x_{\mathcal{R}}$','Interpreter','latex')
xlabel(y1_label)
ylabel(y2_label)
zlabel(flux_label)
view(ax_xRemax,0,90)
%% Fmax surface
fig_Fmax = figure('Name','Max Fidelity');
clf(fig_Fmax)
ax_Fmax = axes(fig_Fmax);

fsurf(ax_Fmax, Fmax, [0 1], 'EdgeColor','none')
colorbar('FontName', 'Times New Roman')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

zlim([0 1]);
xlabel(y1_label)
ylabel(y2_label)
title('$\mathcal{F}_{max}=\frac{1}{4}\left(1+\frac{3}{4\sqrt{y_1y_2}+2(y_1+y_2)+1}\right)$','Interpreter', 'latex')
view(ax_Fmax,0,90);
%% x_Fmax surface
fig_xFmax = figure('Name','X of Max Fidelity');
clf(fig_xFmax)
ax_xFmax = axes(fig_xFmax);
fsurf(ax_xFmax, x_Fmax, [0 1],'EdgeColor','none')

colorbar('FontName','Times New Roman')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

title("$x_\mathcal{F}=2\sqrt{y_1y_2}$",'Interpreter','latex')
xlabel(y1_label)
ylabel(y2_label)
zlabel(flux_label)
view(0,90);
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

xlabel(y1_label)
ylabel(y2_label)
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

fig_tile = figure();
t = tiledlayout(fig_tile,2,2);

ax2_F = subplot(2,2,1,'parent',fig_tile);
ax2_RE = subplot(2,2,3,'parent',fig_tile);

F_h = matlabFunction(F,'Vars',{x, yA, yB});

Y1 = linspace(0,1,11);
Y2 = linspace(0,1,11);

YS = [0 0;
      0 0.3;
      0.1  0.3;
      0.24 0.24];
domain = [0 2];

hold(ax_F,'on')
hold(ax_RE,'on')
for i = 1:size(YS,1)
    ya = YS(i,1);
    yb = YS(i,2);
    
    fplot(ax_F, @(x) F_h(x,ya,yb), domain, 'DisplayName', sprintf("$(y_1,y_2)=(%.2f,%.2f)$",ya,yb)...
        ,'LineWidth',1.5)
    fplot(ax_RE,@(x) RE(x,ya,yb), domain, 'DisplayName', sprintf("$(y_1,y_2)=(%.2f,%.2f)$",ya,yb)...
        ,'LineWidth',1.5)
end
yl = yline(ax_F,0.5,'--','LineWidth',1.5);
yl.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel(ax_RE,'$x$')
xlabel(ax_F,'$x$')
% ylabel(ax_F,'Fidelity')
% ylabel(ax_RE, 'Dimensionless EBR')
title(ax_F,'Fidelity $\mathcal{F}$')
title(ax_RE,'Dimensionless EBR $\mathcal{R}$')
ylim(ax_RE,[0, 00.65])
legend(ax_F)
legend(ax_RE,'Location','bestoutside')
set(findall(fig_F,'-property','FontSize'),'FontSize',18);
set(findall(fig_RE,'-property','FontSize'),'FontSize',18);

% Copying to subplots
axcp_F = copyobj(ax_F,fig_tile);
set(axcp_F, 'Position',get(ax2_F,'position'));
delete(ax2_F);

axcp_RE = copyobj(ax_RE,fig_tile);
set(axcp_RE,'Position',get(ax2_RE,'position'));
delete(ax2_RE);

legend(axcp_RE)

set(findall(fig_tile,'-property','FontSize'),'FontSize',18);
% for i = 1:numel(Y2)
%     yb = Y2(i);
%     b_gifoverwrite = i==1;
%     for ya = Y1
%         fplot(ax_F, @(x)F_h(x, ya, yb), domain, 'DisplayName', sprintf("y_1=%.2f",ya));
%         title(ax_F,sprintf("y_2=%.2f",yb))
%         
%         fplot(ax_RE, @(x)RE(x, ya, yb), domain, 'DisplayName', sprintf("y_1=%.2f",ya));
%         title(ax_RE,sprintf("y_2=%.2f",yb))
%         
%         hold(ax_F,'on')
%         hold(ax_RE,'on')
%     end
%     
%     ylim(ax_F, [0.25 1])
%     ylim(ax_RE, [0 1])
%     
%     yticks(ax_F,[0.25 0.3:0.1:1])
%     xlabel(ax_F,flux_label)
%     xlabel(ax_RE,flux_label)
%     
%     ylabel(ax_F,'Fidelity')
%     ylabel(ax_RE,'Normalized EBR')
%     
%     legend(ax_F)
%     legend(ax_RE)
%     
%     hold(ax_F,'off')
%     hold(ax_RE,'off')
%     
%     writeToGIF(fig_F,gif_F,b_gifoverwrite, 'DelayTime',0.75)
%     writeToGIF(fig_RE,gif_RE,b_gifoverwrite,'DelayTime',0.75)
% end

% close(fig_F)
% close(fig_RE)
%% Plot boundary region
u = linspace(0,1,100);
v = y_nzRE(u,D,0);
fig_bounds = figure();
clf(fig_bounds);
ax_bounds = axes(fig_bounds);

color = 'r';
area(u,v,'FaceColor',color,'FaceAlpha',0.2,'EdgeColor',color,'LineWidth',1.5,'LineStyle','--');
title("Region where $\mathcal{R}_{max}>0$",'Interpreter','latex')
xlabel("$y_1$")
ylabel("$y_2$")
text(0.3,0.5,"$(y_1-y_2)^2-2(y_1+y_2)+1>0$",'Color',color,'Interpreter','latex')
set(findall(fig_bounds,'-property','FontSize'),'FontSize', 18)
%% Plot boundary contours
XX = REmax./max(REmax,[],'all');
XX = Fmax;
Cs = [-2.5:0.5:0.5, 0:0.25:1];

vv = y_nzRE(u(:),D,Cs(:).');
Z = (u-u.').^2-2.*(u+u.')+1;
fig_boundc = figure();
ax_boundc = axes(fig_boundc);
[C,h] = contour(ax_boundc, u, u.',Z,'LevelList',Cs,'LineColor','k');
t = clabel(C);

hold(ax_boundc,'on')

% Add contour labels and curves to transform group
g = hgtransform('Parent',ax_boundc);
set(h,'Parent',g)
set(t,'Parent',g)

% Shifting contour plots up by 2 to be above surface
M = makehgtform('translate',[0 0 2]);
set(g,'Matrix',M)

% plotting xREmax surface
contour(ax_boundc,Y,Y,XX);%,'EdgeAlpha',0)
colorbar(ax_boundc)
caxis([0,max(XX,[],'all')])

% Modifying contour labels
txts = findall(t,'-property','String');
new_txts = compose('$z=%g$',string(get(txts,'String')));
set(txts,{'String'},new_txts);

view(0,90)

set(findall(fig_boundc,'-property','FontSize'),'FontSize',12);
title(ax_boundc,"\textbf{\boldmath${\mathcal{R}_{max}}$ with contours of \boldmath$z=(y_1-y_2)^2-2(y_1+y_2)+1$}")
xlabel(ax_boundc,'$y_1$')
ylabel(ax_boundc, '$y_2$')
%% Plot Q vs Rmax
fig_QR = figure();
ax_QR = axes(fig_QR);

Q = (YA - YB).^2-2.*(YA+YB)+1;

ind_lin = find(Q >= 0);
Q_lin = Q(ind_lin);
RE_lin = REmax(ind_lin);
sh = scatter(ax_QR,Q_lin,RE_lin,3);
sh.Annotation.LegendInformation.IconDisplayStyle = 'off';
title('\textbf{Scatter Plot of \boldmath$\mathcal{R}_{max}$ vs. \boldmath$Q=(y_1-y_2)^2-2(y_1+y_2)+1$}')
xlabel('$Q$')
ylabel('$\mathcal{R}_{max}$')


% Perform fit
fit_degree = 1;
coeff = polyfit(Q_lin, RE_lin, fit_degree);

% Generate polynomial fit string
coeff_format_spec = arrayfun(@(n)['%gx^',num2str(n)],1:fit_degree, 'UniformOutput', false);
coeff_format_spec = [coeff_format_spec(:) ; '%g'];
coeff_str =compose(strjoin(coeff_format_spec,'+'),coeff);
coeff_str = ['$',coeff_str{:},'$'];

% Generate polynomial function handle
fit = @(x) sum(coeff(:).*x(:).');
hold(ax_QR,'on')
fplot(fit,[0 1],'DisplayName', string(coeff_str),'LineWidth',1.5,'LineStyle',':');

legend
