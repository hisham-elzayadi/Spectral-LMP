function plot_eigenvector_diagnostics(MixingMatrix,LeakageMatrix,tauActual,eta,delta,k,n)

%% Plot heat maps

figure;
clf;

tLeak = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Fixed logarithmic color scale
logMin = -16;
logMax = 0;

% Avoid log10(0)
plotFloor = 10^logMin;

MixingLog = log10(max(abs(MixingMatrix),plotFloor));
LeakageLog = log10(max(abs(LeakageMatrix),plotFloor));

% Use the same colormap for both heat maps
map = parula(256);

%% Heat map 1: dominant-eigenvector correspondence

ax1 = nexttile(tLeak,1);

imagesc(ax1,MixingLog);

axis(ax1,'tight');
axis(ax1,'xy');
axis(ax1,'square');

colormap(ax1,map);
clim(ax1,[logMin logMax]);

ax1.FontSize = 20;
ax1.TickLabelInterpreter = 'latex';
ax1.LineWidth = 1;

box(ax1,'on');

xlabel(ax1,'$v_i$', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel(ax1,'$\widehat v_i$', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Rotation',0);

title(ax1,'$\left|\widehat V_k^TQ_k\right|$', ...
    'Interpreter','latex', ...
    'FontSize',20);

% Color bar
cb1 = colorbar(ax1);
cb1.TickLabelInterpreter = 'latex';
cb1.FontSize = 20;

cbTicks = logMin:2:logMax;
cb1.Ticks = cbTicks;

cb1.TickLabels = arrayfun( ...
    @(p) sprintf('$10^{%d}$',p), ...
    cbTicks, ...
    'UniformOutput',false);

%% Heat map 2: complement-subspace contamination

ax2 = nexttile(tLeak,2);

imagesc(ax2,LeakageLog);

axis(ax2,'tight');
axis(ax2,'xy');
axis(ax2,'square');

colormap(ax2,map);
clim(ax2,[logMin logMax]);

ax2.FontSize = 20;
ax2.TickLabelInterpreter = 'latex';
ax2.LineWidth = 1;

box(ax2,'on');

xlabel(ax2,'$v_i$', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel(ax2,'$\widehat v_i$', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Rotation',0);

title(ax2,'$\left|\widehat V_k^TQ_{k+1:n}\right|$', ...
    'Interpreter','latex', ...
    'FontSize',20);

% Color bar
cb2 = colorbar(ax2);
cb2.TickLabelInterpreter = 'latex';
cb2.FontSize = 20;

cb2.Ticks = cbTicks;

cb2.TickLabels = arrayfun( ...
    @(p) sprintf('$10^{%d}$',p), ...
    cbTicks, ...
    'UniformOutput',false);

% Show original complementary eigenvector indices k+1,...,n
numberOfTicks = 4;

tickPositions = round( ...
    linspace(1,n-k,numberOfTicks));

tickLabels = k + tickPositions;

xticks(ax2,tickPositions);
xticklabels(ax2,string(tickLabels));

% Figure size
set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 28 21]);

% Optional export
% exportgraphics(gcf,'ritzit_Diagnostics.pdf', ...
%     'ContentType','vector');

%% Plot actual eigenvector error and its components

figure;
clf;

ax3 = axes;

% Total eigenvector error
hTau = semilogy(ax3,1:k,tauActual, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',10, ...
    'DisplayName', ...
    '$\tau_i=\| \Delta v_i\|_2$');

hold(ax3,'on');

% Complement-subspace component
hEta = semilogy(ax3,1:k,eta, ...
    '--s', ...
    'LineWidth',2, ...
    'MarkerSize',5, ...
    'DisplayName', ...
    '$\eta_i=\|Q_{k+1:n}^T \Delta v_i\|_2$');

% Dominant-subspace component
hDelta = semilogy(ax3,1:k,delta, ...
    '-.^', ...
    'LineWidth',2, ...
    'MarkerSize',5, ...
    'DisplayName', ...
    '$\delta_i=\|Q_k^T \Delta v_i\|_2$');

hold(ax3,'off');

%% Axes formatting

grid(ax3,'on');
box(ax3,'on');

xlabel(ax3,'$\widehat v_i$', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel(ax3,'Error', ...
    'Interpreter','latex', ...
    'FontSize',20);

legend(ax3,[hTau,hEta,hDelta], ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',20, ...
    'Box','off');

ax3.FontSize = 20;
ax3.TickLabelInterpreter = 'latex';
ax3.LineWidth = 1;

xlim(ax3,[1 k]);

%% Automatically determine logarithmic y-axis limits

yl = ylim(ax3);

emin = floor(log10(yl(1)));
emax = ceil(log10(yl(2)));

% Guarantee at least three major logarithmic ticks
if emax - emin < 2
    mid = round((emin + emax)/2);
    emin = mid - 1;
    emax = mid + 1;
end

ylim(ax3,[10^emin 10^emax]);
yticks(ax3,10.^(emin:emax));

ax3.XGrid = 'on';
ax3.YGrid = 'on';
ax3.YMinorGrid = 'off';

% Figure size
set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 24 16]);

% Optional export
% exportgraphics(gcf,'ritzit_IndividualErrors.pdf', ...
%     'ContentType','vector');

end