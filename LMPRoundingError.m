%% Spectral-LMP rounding-error experiments: dominant and complement subspaces
%
% This experiment investigates the effect of the cluster point theta on
% rounding errors in the application of the Spectral-LMP. The observed
% relative errors are compared with the corresponding theoretical bounds
% for vectors in the dominant and complement subspaces. The experiment
% also illustrates the cluster points that minimize the derived bounds.
    
%% Problem setup

% Fix the random seed for reproducibility
rng(1,'twister');

% Generate an n-by-n orthogonal matrix Q whose columns represent eigenvectors
n = 5000;
A = 10*randn(n);
[Q,~] = qr(A,"econ");

% Number of dominant eigenpairs used in the Spectral-LMP
k = 100;

%% Generate eigenvalues and cluster points for the dominant subspace

% Dominant eigenvalues lambda_1,...,lambda_k
lamdas = logspace(16,13,100);

% Harmonic-mean cluster point for the dominant-subspace bound
w=1./lamdas;
W=sum(w);
cumweightes=cumsum(w);
weightedmedian=find(cumweightes >= W/2,1);
thetaWM=lamdas(weightedmedian);
% Candidate values of theta spanning a wide range
eigens = [logspace(19,17,k) lamdas logspace(12,0,k)];

%% Compute the bound and observed error in the dominant subspace

% Unit roundoff and gamma_n
u = eps;
gamma = n*u/(1-n*u);

% Preallocate theoretical bound and observed relative error
bound = zeros(1,3*k);
e = zeros(1,3*k);

% Sort candidate cluster points
eigens = sort(eigens,'ascend');

% Locate the cluster points/eigenvalues highlighted in the plot
thetaWMIndex = find(eigens == thetaWM);
lamda_k = find(eigens == lamdas(k));
lamda_n = find(eigens == 1);
lamda_1 = find(eigens == 1e16);

% Evaluate the theoretical dominant-subspace rounding-error bound
% for each candidate cluster point theta
for i = 1:length(eigens)
  
    theta = eigens(i);
    % alpha_j = 1 - theta/lambda_j
    alpha = 1 - theta./lamdas(1:k);

    s = k*u + k*gamma*sum(abs(alpha));

    if theta >= lamdas(1)

        % theta >= lambda_1
        bound(i) = (theta/lamdas(k))*s;

    elseif theta >= lamdas(k)

        % lambda_k <= theta < lambda_1
        bound(i) = (lamdas(1)/lamdas(k))*s;

    else

        % 0 < theta < lambda_k
        bound(i) = (lamdas(1)/theta)*s;

    end
end
% Generate random coefficients and construct a vector lying entirely
% in the dominant eigenspace span{v_1,...,v_k}
q = normrnd(0,1e3/3,1,k);
x = Q(:,1:k)*q';

% Compare the computed Spectral-LMP action with the exact
% dominant-subspace action
for i = 1:3*k

    % Computed application of the Spectral-LMP
    y = Spectral_LMP(Q,lamdas,k,eigens(i),x);

    % Exact action in the dominant subspace:
    % H_k y_D = \sum (theta*\beta_i/lambda_i) v_i
    ytrue = Q(:,1:k)*((eigens(i)*q./lamdas(1:k)))';

    % Observed relative forward error
    e(i) = norm(ytrue-y)/ norm(eigens(i)*(q./lamdas(1:k)));
end

%% Generate the plot for the dominant subspace

figure(1);
clf;
c = colororder;

t = tiledlayout(2,1,...
    'TileSpacing','compact',...
    'Padding','compact');

% Common y-axis label for the two panels
ylabel(t,...
    '$\displaystyle\frac{\|H_ky_D-\mathrm{fl}(H_ky_D)\|_2}{\|H_ky_D\|_2}$',...
    'Interpreter','latex',...
    'FontSize',20);

% Common x-axis label
xlabel(t,...
    '$\theta$',...
    'Interpreter','latex',...
    'FontSize',20);

% ==================== (a) Bound ====================

ax1 = nexttile;
hold(ax1,'on');

% Theoretical rounding-error bound
hBound = plot(ax1,eigens,bound,...
    'r-','LineWidth',2.5);

% Dummy blue curve used only to include the observed error
% in the shared legend
hObserved = plot(ax1,nan,nan,...
    'b-','LineWidth',2.5);

% Harmonic-mean cluster point
hWM = plot(ax1,thetaWM,bound(thetaWMIndex),...
    'o','color',c(3,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_1
hL1 = plot(ax1,eigens(lamda_1),bound(lamda_1),...
    'd','color',c(4,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_k
hLk = plot(ax1,eigens(lamda_k),bound(lamda_k),...
    'd','color',c(5,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_n
hLn = plot(ax1,eigens(lamda_n),bound(lamda_n),...
    'p','color',c(6,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

hold(ax1,'off');

% Use logarithmic scales on both axes
ax1.YScale = 'log';
ax1.XScale = 'log';

grid(ax1,'on');
box(ax1,'on');

ax1.FontSize = 20;
ax1.LineWidth = 1;
ax1.TickLabelInterpreter = 'latex';

% Suppress x-axis tick labels in the upper panel
ax1.XTickLabel = [];

title(ax1,'(a) Theoretical bound',...
    'Interpreter','latex',...
    'FontSize',20,...
    'HorizontalAlignment','left');

% Shared legend for the bound and observed-error panels
lgd1 = legend(ax1,...
    [hBound,hObserved,hWM,hL1,hLk,hLn],...
    {'Bound',...
     'Observed error',...
     '$\theta_{WMr}$',...
     '$\lambda_1$',...
     '$\lambda_k$',...
     '$\lambda_n$'},...
    'Interpreter','latex',...
    'Orientation','horizontal',...
    'NumColumns',3,...
    'FontSize',20,...
    'Box','off',...
    'Location','northoutside');

% ==================== (b) Observed error ====================

ax2 = nexttile;
hold(ax2,'on');

% Observed relative forward error
plot(ax2,eigens,e,...
    'b-','LineWidth',2.5);

% Harmonic-mean cluster point
plot(ax2,thetaWM,e(thetaWMIndex),...
    'o','color',c(3,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_1
plot(ax2,eigens(lamda_1),e(lamda_1),...
    'd','color',c(4,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_k
plot(ax2,eigens(lamda_k),e(lamda_k),...
    'd','color',c(5,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_n
plot(ax2,eigens(lamda_n),e(lamda_n),...
    'p','color',c(6,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

hold(ax2,'off');

% Use logarithmic scales on both axes
ax2.YScale = 'log';
ax2.XScale = 'log';

grid(ax2,'on');
box(ax2,'on');

ax2.FontSize = 20;
ax2.LineWidth = 1;
ax2.TickLabelInterpreter = 'latex';

title(ax2,'(b) Observed error',...
    'Interpreter','latex',...
    'FontSize',20,...
    'HorizontalAlignment','left');

% Use the same y-axis ticks in both panels
ticks = 10.^(-15:5:5);
yticks(ax1,ticks);
yticks(ax2,ticks);

% Set figure dimensions
set(gcf,...
    'Color','w',...
    'Units','centimeters',...
    'Position',[2 2 24 16]);

% Export the dominant-subspace figure as a vector PDF
exportgraphics(gcf,...
    'rounding_error_Dom.pdf',...
    'ContentType','vector',...
    'BackgroundColor','white');

%% ==================== Complement subspace ====================

% Cluster point minimizing the complement-subspace bound
thetaC = sum(1./lamdas(1:k))/ sum(1./lamdas(1:k).^2);

% Candidate values of theta spanning a wide range
eigens = [logspace(19,17,k) lamdas logspace(12,0,k) thetaC];

%% Compute the bound and observed error in the complement subspace

% Sort candidate cluster points
eigens = sort(eigens,'ascend');

% Locate the cluster point and eigenvalues highlighted in the plot
thetaCIndex = find(eigens == thetaC);
lamda_k = find(eigens == lamdas(k));
lamda_n = find(eigens == 1);
lamda_1 = find(eigens == 1e16);

% Evaluate the theoretical complement-subspace rounding-error bound
for i = 1:3*k+1

    % Spectral-LMP coefficients alpha_i = 1 - theta/lambda_i
    alpha = 1 - (eigens(i)./lamdas(1:k));

    % Complement-subspace bound
    bound(i) = gamma*norm(alpha,2) + k*u;
end

% Generate random coefficients and construct a vector lying entirely
% in the complement of the dominant eigenspace
q = normrnd(0,1e3/3,1,n-k);
x = Q(:,k+1:n)*q';

% In exact arithmetic, H_k acts as the identity on the complement.
% Measure the deviation of the computed result from x.
for i = 1:3*k+1

    % Computed application of the Spectral-LMP
    y = Spectral_LMP(Q,lamdas,k,eigens(i),x);

    % Observed relative forward error
    e(i) = norm(x-y,2)/norm(x,2);
end

%% Generate the plot for the complement subspace

figure(2);
clf;
c = colororder;

p = tiledlayout(2,1,...
    'TileSpacing','compact',...
    'Padding','compact');

% Common y-axis label for the two panels
ylabel(p,...
    '$\displaystyle\frac{\|H_ky_C-\mathrm{fl}(H_ky_C)\|_2}{\|H_ky_C\|_2}$',...
    'Interpreter','latex',...
    'FontSize',20);

% Common x-axis label
xlabel(p,...
    '$\theta$',...
    'Interpreter','latex',...
    'FontSize',20);

% ==================== (a) Bound ====================

ax1 = nexttile;
hold(ax1,'on');

% Theoretical rounding-error bound
hBound = plot(ax1,eigens,bound,...
    'r-','LineWidth',2.5);

% Dummy blue curve used only to include the observed error
% in the shared legend
hObserved = plot(ax1,nan,nan,...
    'b-','LineWidth',2.5);

% Complement-subspace minimizing cluster point
hWAM = plot(ax1,thetaC,bound(thetaCIndex),...
    'o','color',c(3,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_1
hL1 = plot(ax1,eigens(lamda_1),bound(lamda_1),...
    'd','color',c(4,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_k
hLk = plot(ax1,eigens(lamda_k),bound(lamda_k),...
    'd','color',c(5,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_n
hLn = plot(ax1,eigens(lamda_n),bound(lamda_n),...
    'p','color',c(6,:),...
    'LineStyle','none',...
    'MarkerSize',10,...
    'LineWidth',2.5);

hold(ax1,'off');

% Use logarithmic scales on both axes
ax1.YScale = 'log';
ax1.XScale = 'log';

grid(ax1,'on');
box(ax1,'on');

ax1.FontSize = 20;
ax1.LineWidth = 1;
ax1.TickLabelInterpreter = 'latex';

% Suppress x-axis tick labels in the upper panel
ax1.XTickLabel = [];

title(ax1,'(a) Theoretical bound',...
    'Interpreter','latex',...
    'FontSize',20,...
    'HorizontalAlignment','left');

% Shared legend for the bound and observed-error panels
lgd = legend(ax1,...
    [hBound,hObserved,hWAM,hL1,hLk,hLn],...
    {'Bound',...
     'Observed error',...
     '$\theta_{WAM}$',...
     '$\lambda_1$',...
     '$\lambda_k$',...
     '$\lambda_n$'},...
    'Interpreter','latex',...
    'Orientation','horizontal',...
    'NumColumns',3,...
    'FontSize',20,...
    'Box','off',...
    'Location','northoutside');

% ==================== (b) Observed error ====================

ax2 = nexttile;
hold(ax2,'on');

% Observed relative forward error
plot(ax2,eigens,e,...
    'b-','LineWidth',2.5);

% Complement-subspace minimizing cluster point
plot(ax2,thetaC,e(thetaCIndex),...
    'o','color',c(3,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_1
plot(ax2,eigens(lamda_1),e(lamda_1),...
    'd','color',c(4,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_k
plot(ax2,eigens(lamda_k),e(lamda_k),...
    'd','color',c(5,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

% Marker for lambda_n
plot(ax2,eigens(lamda_n),e(lamda_n),...
    'p','color',c(6,:),...
    'MarkerSize',10,...
    'LineWidth',2.5);

hold(ax2,'off');

% Use logarithmic scales on both axes
ax2.YScale = 'log';
ax2.XScale = 'log';

grid(ax2,'on');
box(ax2,'on');

ax2.FontSize = 20;
ax2.LineWidth = 1;
ax2.TickLabelInterpreter = 'latex';

title(ax2,'(b) Observed error',...
    'Interpreter','latex',...
    'FontSize',20,...
    'HorizontalAlignment','left');

% Use the same y-axis ticks in both panels
ticks = 10.^(-15:5:5);
yticks(ax1,ticks);
yticks(ax2,ticks);

% Set figure dimensions
set(gcf,...
    'Color','w',...
    'Units','centimeters',...
    'Position',[2 2 24 16]);

% Export the complement-subspace figure as a vector PDF
exportgraphics(gcf,...
    'rounding_error_Comp.pdf',...
    'ContentType','vector',...
    'BackgroundColor','white');