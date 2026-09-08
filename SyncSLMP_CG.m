%% Synthetic Spectral-LMP experiments
%
% Three synthetic perturbation experiments are performed. The experiments
% differ only in the prescribed perturbation levels tau and the parameter
% gamma controlling the complement-subspace component of the perturbation.

%% Common problem setup

rng(1,'twister');

n = 5000;

A = 10*randn(n);
[Q,~] = qr(A,"econ");

k = 100;

lambdas = [linspace(1e16,1e14,k) linspace(1e12,1,n-k)];

Afun = @(z) Q*(lambdas'.*(Q'*z));

x = randn(n,1);
b = Afun(x);

%% Define the synthetic experiments

experiments(1).tau   = linspace(1e-6,1e0,k);
experiments(1).gamma = 0;
experiments(1).name  = 'Synthetic 1';

experiments(2).tau   = linspace(1e0,1e-16,k);
experiments(2).gamma = 1;
experiments(2).name  = 'Synthetic 2';

experiments(3).tau   = linspace(1e-12,1e-6,k);
experiments(3).gamma = 0;
experiments(3).name  = 'Synthetic 3';

%% Run the experiments

for experimentID = 1:length(experiments)

    tau   = experiments(experimentID).tau;
    gamma = experiments(experimentID).gamma;

    fprintf('\n========================================\n');
    fprintf('%s\n',experiments(experimentID).name);
    fprintf('gamma = %.2e\n',gamma);
    fprintf('========================================\n');

    %% Generate perturbed eigenvectors

    % Use the same random perturbation matrix in every experiment
    rng(2,'twister');

    E = randn(n,k);

    % Control the component of the perturbation lying in the
    % complementary eigenspace
    E = E - gamma*Q(:,k+1:n)*Q(:,k+1:n)'*E;

    % Normalize the perturbation directions
    E = E./vecnorm(E);

    % Apply the prescribed perturbation levels
    V_hat = Q(:,1:k) + E.*tau;

    % Reorthogonalize the perturbed eigenvectors
    [V_hat,~] = qr(V_hat,0);
    %% Align the signs of approximate and exact eigenvectors
    % 
    for i = 1:k
        if V_hat(:,i)'*Q(:,i) < 0
            V_hat(:,i) = -V_hat(:,i);
        end
    end

    % Actual eigenvector errors after reorthogonalization
    tauActual = vecnorm(V_hat-Q(:,1:k));

    %% Run CG/Spectral-LMP experiment

    results = run_LMP_experiment(Afun,b,x,Q,V_hat,lambdas,[],k,tauActual);

    %% Plot convergence

    plot_convergence(results);

    %% Eigenvector mixing and complement-subspace contamination

    Qk    = Q(:,1:k);
    Qcomp = Q(:,k+1:n);

    % Mixing within the dominant eigenspace
    MixingMatrix = V_hat'*Qk;

    % Leakage into the complementary eigenspace
    LeakageMatrix = V_hat'*Qcomp;

    % Rowwise complementary-subspace leakage
    eta = vecnorm(LeakageMatrix,2,2);

    % Dominant-subspace component of the eigenvector error
    delta = sqrt(max(tauActual'.^2-eta.^2,0));

    % Global complement-subspace contamination
    leakageNorm = norm(LeakageMatrix,2);

    % Largest principal angle
    thetaMax    = asin(min(1,leakageNorm));
    thetaMaxDeg = rad2deg(thetaMax);

    fprintf('\nComplementary-subspace contamination diagnostics\n');
    fprintf('-------------------------------------------------\n');
    fprintf('||V_hat'' Q_comp||_2    = %.6e\n',leakageNorm);
    fprintf('Largest principal angle = %.6e radians\n',thetaMax);
    fprintf('Largest principal angle = %.6e degrees\n',thetaMaxDeg);

    %% Plot eigenvector diagnostics

    plot_eigenvector_diagnostics(MixingMatrix,LeakageMatrix,tauActual,eta,delta,k,n);

    %% Store experiment results

    experiments(experimentID).results = results;
    experiments(experimentID).tauActual = tauActual;
    experiments(experimentID).eta = eta;
    experiments(experimentID).delta = delta;
    experiments(experimentID).leakageNorm = leakageNorm;
    experiments(experimentID).thetaMax = thetaMax;

end