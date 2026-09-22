%% Exact spectral clustering vs. Spectral-LMP in finite precision
%
% This experiment compares CG applied to an exactly clustered system with
% PCG using the practical Spectral-LMP. Both approaches use the same exact
% spectral information and produce the same clustered spectrum in exact
% arithmetic. Therefore, differences in their convergence behavior reveal
% the effect of rounding errors introduced when applying the Spectral-LMP.

%% Generate an orthogonal eigenvector matrix

% Fix the random seed for reproducibility
rng(1,'twister');

n = 5000;

% Generate a random matrix and compute an orthogonal basis Q
A = 10*randn(n);
[Q,~] = qr(A,"econ");

% Number of dominant eigenvalues modified by the Spectral-LMP
k = 100;

%% Prepare the linear systems

% Construct the spectrum of the original SPD matrix.
% The first k eigenvalues form the dominant part of the spectrum.
lamdas = [linspace(1e8,1e6,100) linspace(1e5,1,4900)];

% Construct the exactly clustered spectrum by replacing the first k
% eigenvalues with the cluster point theta = 1.
lamdasPre = [linspace(1,1,100) linspace(1e5,1,4900)];

% Matrix-vector product with the original matrix
% A = Q*diag(lamdas)*Q'
A = @(z) Q*(lamdas'.*(Q'*z));

% Matrix-vector product with the exactly clustered matrix
% APre = Q*diag(lamdasPre)*Q'
APre = @(z) Q*(lamdasPre'.*(Q'*z));

% Generate an exact solution
x = randn(n,1);

% Construct the corresponding right-hand side for the original system
b = A(x);

% Construct the right-hand side for the exactly clustered system so that
% both systems have the same exact solution x
bPre = APre(x);

%% CG and PCG experiments

% Apply standard CG to the original unpreconditioned system
[~,~,~,~,~,XCG] = pcg_modified(A,b,1e-16,30);

% Define the Spectral-LMP using the exact dominant eigenvectors and
% eigenvalues, with cluster point theta = 1
LMP = @(u) Spectral_LMP(Q,lamdas,k,1,u);

% Apply PCG to the original system using the practical Spectral-LMP
[~,~,~,~,~,XLMP] = pcg_modified(A,b,1e-16,30,LMP);

% Apply CG directly to the exactly clustered system.
% In exact arithmetic, this system is spectrally equivalent to the
% Spectral-LMP preconditioned system.
[~,~,~,~,~,X] = pcg_modified(APre,bPre,1e-16,30);

%% Compute relative solution errors

% Relative error history for standard CG
eCG = XCG - x;
xCGerr2 = vecnorm(eCG)/norm(x);

% Relative error history for PCG with Spectral-LMP
e1 = XLMP - x;
xLMPerr2 = vecnorm(e1)/norm(x);

% Relative error history for CG applied to the exactly clustered system
e = X - x;
xerr2 = vecnorm(e)/norm(x);

%% Plot the convergence histories

figure(2);
clf;

ax = axes;
hold(ax,'on');

% Iteration indices for the three methods
iterCG    = 0:size(xCGerr2,2)-1;
iterExact = 0:size(xerr2,2)-1;
iterLMP   = 0:size(xLMPerr2,2)-1;

% Standard CG on the original system
hCG = semilogy(ax,iterCG,xCGerr2,...
    '-','LineWidth',2.5);

% CG on the exactly clustered system
hExact = semilogy(ax,iterExact,xerr2,...
    '-d','LineWidth',2.5);

% PCG using the practical Spectral-LMP
hLMP = semilogy(ax,iterLMP,xLMPerr2,...
    '--*',...
    'LineWidth',2.5,...
    'MarkerSize',10,...
    'MarkerIndices',1:2:length(iterLMP));

hold(ax,'off');

%% Axes formatting

grid(ax,'on');
box(ax,'on');

ax.FontSize = 20;
ax.LineWidth = 1;
ax.TickLabelInterpreter = 'latex';
ax.YScale = 'log';

% Leave one extra iteration of space after the longest curve
lastIter = max([iterExact,iterLMP]);
xlim(ax,[0 lastIter+1]);

% Use integer iteration ticks
xticks(ax,0:2:lastIter+1);

xlabel(ax,'CG iterations',...
    'Interpreter','latex',...
    'FontSize',20);

% Relative forward error with respect to the exact solution
ylabel(ax,...
    '$\displaystyle\frac{\|x^\ast-\hat{x}_\ell\|_2}{\|x^\ast\|_2}$',...
    'Interpreter','latex',...
    'FontSize',20,...
    'Rotation',90);

ylim(ax,[1e-1 2e1]);

%% Legend

lgd = legend(ax,...
    [hCG,hExact,hLMP],...
    {'CG','Exact spectral clustering',...
     'Spectral-LMP'},...
    'Interpreter','latex',...
    'Orientation','horizontal',...
    'FontSize',20,...
    'Box','off',...
    'Location','northoutside');

%% Figure dimensions and export

set(gcf,...
    'Color','w',...
    'Units','centimeters',...
    'Position',[2 2 24 16]);
