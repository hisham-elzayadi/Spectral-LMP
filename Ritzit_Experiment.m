%% Ritzit experiment
%
% This experiment studies the performance of the Spectral-LMP when the
% dominant eigenvectors are estimated using the ritzit
% method. The convergence of CG with different cluster points is compared,
% and the accuracy and complementary-subspace contamination of the
% approximate dominant eigenvectors are examined.

%% Prepare the eigenvectors

% Fix the random seed for reproducibility
rng(1,'twister');

n = 5000;

% Generate an orthogonal eigenvector matrix Q
A = 10*randn(n);
[Q,~] = qr(A,"econ");

%% Prepare the linear system

% Number of dominant eigenpairs
k = 100;

% Prescribed spectrum of the SPD matrix
lambdas = [linspace(1e16,1e14,k) linspace(1e12,1,n-k)];

% Matrix-vector product with
% A = Q*diag(lambdas)*Q'
Afun = @(z) Q*(lambdas'.*(Q'*z));

% Generate the exact solution and corresponding right-hand side
x = randn(n,1);
b = Afun(x);

%% Estimate the dominant eigenpairs using randomized Nyström

% Fix the random seed used by the randomized eigensolver
rng(2,'twister');

% Compute approximate dominant eigenvectors and eigenvalues
[V_hat,Theta] = REVD_ritzit(Afun,k,5,n);

% Extract the approximate dominant eigenvalues.
% Use diag(Theta) here if Theta is returned as a diagonal matrix.
lambdas_D = diag(Theta);

%% Align approximate and exact eigenvectors

% Resolve the sign ambiguity of the approximate eigenvectors so that
% each approximate eigenvector has the same orientation as its exact
% counterpart.
for i = 1:k
    if V_hat(:,i)'*Q(:,i) < 0
        V_hat(:,i) = -V_hat(:,i);
    end
end

% Actual error in each approximate dominant eigenvector
tauActual = vecnorm(V_hat-Q(:,1:k));

%% Run CG/Spectral-LMP experiment

% Run CG and Spectral-LMP using the approximate Nyström eigenvectors
% and approximate dominant eigenvalues
results = run_LMP_experiment(Afun,b,x,Q,V_hat,lambdas,lambdas_D,k,tauActual);

%% Plot convergence results

plot_convergence(results);

%% Eigenvector mixing and complement-subspace contamination

% Exact dominant and complementary eigenspaces
Qk    = Q(:,1:k);
Qcomp = Q(:,k+1:n);

% Mixing within the dominant eigenspace:
% M(i,j) = v_hat_i^T v_j,   i,j = 1,...,k
MixingMatrix = V_hat'*Qk;

% Leakage into the complementary eigenspace:
% L(i,j) = v_hat_i^T v_{k+j}
LeakageMatrix = V_hat'*Qcomp;

% Rowwise complementary-subspace leakage:
% eta_i = ||Qcomp^T v_hat_i||_2
eta = vecnorm(LeakageMatrix,2,2);

% Dominant-subspace component of the eigenvector error
delta = sqrt(max(tauActual'.^2-eta.^2,0));

% Global complementary-subspace contamination:
% ||V_hat^T Qcomp||_2 = sin(theta_max)
leakageNorm = norm(LeakageMatrix,2);

% Largest principal angle between the exact and approximate
% dominant eigenspaces
thetaMax = asin(min(1,leakageNorm));
thetaMaxDeg = rad2deg(thetaMax);

%% Display eigenvector diagnostics

fprintf('\nNyström eigenvector diagnostics\n');
fprintf('---------------------------------------------\n');
fprintf('||V_hat'' Q_comp||_2    = %.6e\n',leakageNorm);
fprintf('Largest principal angle = %.6e radians\n',thetaMax);
fprintf('Largest principal angle = %.6e degrees\n',thetaMaxDeg);

%% Plot eigenvector diagnostics

plot_eigenvector_diagnostics(MixingMatrix,LeakageMatrix,tauActual,eta,delta,k,n);