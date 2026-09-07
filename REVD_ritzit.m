function [U, Theta] = REVD_ritzit(A, k, ell, n)
% REVD_RITZIT Randomised eigenvalue decomposition based on ritzit
%
% Inputs:
%   A   : symmetric matrix or function handle Y = A(X)
%   k   : target rank
%   ell : oversampling parameter
%   n   : dimension, required only when A is a function handle
%
% Outputs:
%   U     : approximate eigenvectors, size n-by-k
%   Theta : approximate largest eigenvalues, k-by-k diagonal matrix

    if isnumeric(A)
        n = size(A,1);
        applyA = @(X) A * X;
    else
        if nargin < 4
            error('For a function handle, specify n.');
        end
        applyA = @(X) A(X);
    end

    m = k + ell;

    % 1. Gaussian random matrix
    G = randn(n, m);

    % 2. Orthonormalize G
    [G3, ~] = qr(G, 0);

    % 3. Sample matrix
    Y3 = applyA(G3);

    % 4. QR decomposition
    [Z3, R3] = qr(Y3, 0);

    % 5. Form K3 = R3 R3^T
    K3 = R3 * R3';

    % 6. Eigenvalue decomposition
    [W3, D] = eig(K3);

    % Sort eigenvalues in decreasing order
    [eigvals_sq, idx] = sort(diag(D), 'descend');
    W3 = W3(:, idx);

    % Approximate eigenvalues of A
    eigvals = sqrt(max(eigvals_sq, 0));

    % 7. Keep first k eigenvalues
    eigvals = eigvals(1:k);

    % 8. Keep first k Ritz vectors
    W3 = W3(:, 1:k);

    % 9. Form approximate eigenvectors
    U = Z3 * W3;

    % Diagonal matrix of approximate eigenvalues
    Theta = diag(eigvals);
end