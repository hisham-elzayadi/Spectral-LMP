function results = run_LMP_experiment(Afun,b,x,Q,V,lambdas,lambdas_D,k,tauActual)

if isempty(lambdas_D)
    lambdas_D = lambdas(1:k);
end

tol = 1e-16;
Maxit = 50;

%% Unpreconditioned CG

[~,~,~,~,~,X] = pcg_modified(Afun,b,tol,Maxit);

%% Weighted-median cluster point

w = 2*tauActual + tauActual.^2;
W = sum(w);
cumweights = cumsum(w);
weightedmedian = find(cumweights >= W/2,1);

thetaWMs = lambdas_D(weightedmedian);

LMPWMs = @(u) Spectral_LMP(V,lambdas_D,k,thetaWMs,u);

[~,~,~,~,~,XWMs] = pcg_modified(Afun,b,tol,Maxit,[],LMPWMs);

%% rounding error weighted median cluster point
w = 1 ./ lambdas_D;
W=sum(w);
cumweightes=cumsum(w);
weightedmedian=find(cumweightes >= W/2,1);
thetaWMr=lambdas_D(weightedmedian);
LMPWMr= @(u) Spectral_LMP(V,lambdas_D,k,thetaWMr,u);

[~,~,~,~,~,XWMr] = pcg_modified(Afun,b,tol,Maxit,[],LMPWMr);

%% Residual-based cluster point

n = length(lambdas);

eta = lambdas(k+1:end).*(Q(:,k+1:n)'*x)';

den = sum(eta.^2);
num = sum(lambdas(k+1:end).*(eta.^2));

theta_r = num/den;

LMPr = @(u) Spectral_LMP(V,lambdas_D,k,theta_r,u);

[~,~,~,~,~,Xr] = pcg_modified(Afun,b,tol,Maxit,[],LMPr);

%% Midpoint cluster point

theta_m = (lambdas(k+1) + lambdas(n))/2;

LMPm = @(u) Spectral_LMP(V,lambdas_D,k,theta_m,u);

[~,~,~,~,~,Xm] = pcg_modified(Afun,b,tol,Maxit,[],LMPm);

%% theta = lambda_k

theta_k = lambdas_D(k);

LMPk = @(u) Spectral_LMP(V,lambdas_D,k,theta_k,u);

[~,~,~,~,~,Xk] = pcg_modified(Afun,b,tol,Maxit,[],LMPk);

%% theta = lambda_{k+1}

theta_kplus1 = lambdas(k+1);

LMPkplus1 = @(u) Spectral_LMP(V,lambdas_D,k,theta_kplus1,u);

[~,~,~,~,~,Xkplus1] = pcg_modified(Afun,b,tol,Maxit,[],LMPkplus1);

%% Compute relative solution errors

results.CG = vecnorm(X-x)/norm(x);

results.WMs = vecnorm(XWMs-x)/norm(x);

results.WMr = vecnorm(XWMr-x)/norm(x);

results.lambda_k = vecnorm(Xk-x)/norm(x);

results.lambda_kplus1 = vecnorm(Xkplus1-x)/norm(x);

results.theta_r = vecnorm(Xr-x)/norm(x);

results.theta_m = vecnorm(Xm-x)/norm(x);

%% Store cluster points as well

results.theta.WMs = thetaWMs;
results.theta.WMr = thetaWMr;
results.theta.lambda_k = theta_k;
results.theta.lambda_kplus1 = theta_kplus1;
results.theta.r = theta_r;
results.theta.m = theta_m;

results.Maxit = Maxit;

end