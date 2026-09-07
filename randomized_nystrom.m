function [U,Theta] = randomized_nystrom(A,k,l,n)
% Randomized Nyström eigendecomposition for symmetric positive semidefinite A
%
% Inputs:
%   A : symmetric positive semidefinite matrix OR function handle
%   k : target rank
%   l : oversampling parameter
%   n : dimension (only required if A is a function handle)
%
% Outputs:
%   U     : n-by-k approximate eigenvectors
%   Theta : k-by-k diagonal matrix of approximate eigenvalues
%
% Approximation:
%   A ≈ U*Theta*U'

    if isnumeric(A)
        n = size(A,1);
        applyA = @(X) A*X;
    else
        if nargin < 4
            error('For a function handle, specify n.');
        end
        applyA = @(X) A(X);
    end

    m = k + l;

    % Step 1: Gaussian random matrix
    G = randn(n,m);

    % Step 2: sample matrix
    Y = applyA(G);

    % Optional numerical stabilization
    nu = sqrt(n) * eps(norm(Y,'fro'));
    Y = Y + nu*G;

    % Step 3: orthonormalize Y
    [Z,~] = qr(Y,0);

    % Step 4: form E1 = A Z and E2 = Z' A Z
    E1 = applyA(Z);
    E2 = Z'*E1;

    % Symmetrize E2 to remove tiny roundoff asymmetry
    E2 = (E2 + E2')/2;

    % Step 5: Cholesky factorization E2 = C' C
    [C,p] = chol(E2);

    % Extra safeguard if Cholesky fails
    if p ~= 0
        delta = sqrt(m) * eps(norm(E2,'fro'));
        C = chol(E2 + delta*eye(m));
    end

    % Step 6: solve F C = E1
    F = E1 / C;

    % Step 7: SVD of F
    [Ufull,S,~] = svd(F,'econ');

    % Steps 8-9: keep first k components
    U = Ufull(:,1:k);
    Sigma = S(1:k,1:k);

    Theta = Sigma.^2;
end