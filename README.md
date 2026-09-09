# Spectral-LMP

MATLAB implementation and numerical experiments for studying the influence of the cluster point in the **Spectral Limited-Memory Preconditioner (Spectral-LMP)**.

The code accompanies the paper

> H. Elzayyadi and J. M. Tabeart,  
> *The Influence of the Cluster Point on Rounding Errors and Sensitivity in the Spectral Limited-Memory Preconditioner.*

The experiments investigate the finite-precision behavior of Spectral-LMP, its sensitivity to inaccurate spectral information, and the effect of different choices of the cluster point.

## Features

- Comparison of practical Spectral-LMP with exact spectral modification.
- Analysis of rounding errors in the dominant and complementary eigenspaces.
- Synthetic experiments with perturbed dominant eigenvectors.
- Comparison of several choices of the cluster point.
- Approximate spectral information computed using randomized Nyström and REVD-ritzit methods.
- Convergence and eigenspace diagnostics for approximate eigenvectors.

## Installation

Clone the repository:

```bash
git clone https://github.com/hisham-elzayadi/Spectral-LMP.git
cd Spectral-LMP
```

The code is written in MATLAB. Place the repository on the MATLAB path:

```matlab
addpath(pwd)
```

The modified PCG routine `pcg_modified.m` requires the MATLAB helper functions `iterchk` and `iterapp`. These functions are not included in this repository and must be accessible on the MATLAB path.

To check whether MATLAB can locate them, run:

```matlab
which iterchk
which iterapp
```

MATLAB should return the location of each function. If either function cannot be found, make the directory containing the required helper functions available on the MATLAB path, for example:

```matlab
addpath('path\to\helper\functions')
```

You can verify again using `which iterchk` and `which iterapp` before running the experiments.

No compilation or additional package installation is required.

## Getting started

To run a synthetic Spectral-LMP experiment with perturbed dominant eigenvectors, use

```matlab
SyncSLMP_CG
```

This generates the test problem, constructs perturbed spectral information, runs CG/PCG with different cluster-point choices, and produces convergence and eigenvector diagnostics.

Other experiments can be run directly from their corresponding MATLAB scripts as described below.

## Numerical experiments

### Exact spectral modification vs. Spectral-LMP

```matlab
ExactLMPvsLMP
```

Compares unpreconditioned CG, CG applied to an exactly clustered system, and PCG using the practical Spectral-LMP.

The exactly clustered and Spectral-LMP-preconditioned systems have the same clustered spectrum in exact arithmetic. Differences in their observed numerical behavior illustrate finite-precision effects arising from the practical application of the preconditioner.

### Rounding-error experiment

```matlab
LMPRoundingError
```

Investigates rounding errors in the dominant and complementary eigenspaces and compares the observed errors with the corresponding theoretical bounds.

### Synthetic eigenvector perturbations

```matlab
SyncSLMP_CG
```

Runs synthetic experiments with perturbed dominant eigenvectors. The perturbation magnitude and its component in the complementary eigenspace can be controlled independently.

The script runs the three synthetic perturbation experiments defined in `SyncSLMP_CG.m` and produces convergence histories together with eigenvector diagnostics.

### Randomized Nyström experiment

```matlab
Nystrom_Experiment
```

Computes approximate dominant eigenpairs using the randomized Nyström method and studies the resulting Spectral-LMP convergence.

The randomized eigensolver is implemented in `randomized_nystrom.m`.

### REVD-ritzit experiment

```matlab
Ritzit_Experiment
```

Computes approximate dominant eigenpairs using the REVD-ritzit method and performs the corresponding convergence and eigenspace diagnostics.

The randomized eigensolver is implemented in `REVD_ritzit.m`.

## Eigenvector diagnostics

For approximate dominant eigenvectors $\widehat V_k$, the experiments distinguish between the magnitude and direction of the spectral perturbations.

For each approximate eigenvector $\widehat v_i$, the total perturbation is measured by

$$
\tau_i = \|\widehat v_i-v_i\|_2,
$$

with dominant- and complementary-eigenspace components

$$
\eta_i = \|Q_k^T(\widehat v_i-v_i)\|_2,
\qquad
\delta_i = \|Q_{k+1:n}^T(\widehat v_i-v_i)\|_2.
$$

The code also visualizes

$$
|\widehat V_k^TQ_k|
\qquad\text{and}\qquad
|\widehat V_k^TQ_{k+1:n}|,
$$

which indicate mixing within the dominant eigenspace and leakage into the complementary eigenspace, respectively.

## Repository structure

The main routines are:

| File | Description |
|---|---|
| `Spectral_LMP.m` | Applies the Spectral-LMP preconditioner. |
| `run_LMP_experiment.m` | Runs CG/PCG for the different cluster-point choices. |
| `pcg_modified.m` | Modified MATLAB PCG implementation that records the CG iterate history. |
| `plot_convergence.m` | Plots relative forward-error convergence histories. |
| `plot_eigenvector_diagnostics.m` | Produces eigenvector accuracy and eigenspace diagnostics. |
| `randomized_nystrom.m` | Randomized Nyström eigensolver. |
| `REVD_ritzit.m` | REVD-ritzit eigensolver. |

`run_LMP_experiment.m` compares the cluster-point choices

$$
\theta_{WM,e},\qquad
\theta_{WM,1/\lambda},\qquad
\lambda_k,\qquad
\lambda_{k+1},\qquad
\theta_r,\qquad
\theta_m.
$$

## Requirements

The experiments use standard MATLAB dense linear-algebra routines, including

```text
qr
svd
eig
chol
```

`LMPRoundingError.m` currently uses `normrnd`, which requires the **Statistics and Machine Learning Toolbox**. The corresponding normally distributed random vectors can alternatively be generated using `randn`.

The supplied experiments use matrices of dimension $n=5000$. Because a dense $5000\times5000$ orthogonal matrix is generated, some experiments require substantial memory and runtime.

## Reproducibility

The experiment scripts use fixed random-number-generator seeds of the form

```matlab
rng(...,'twister')
```

to make the randomly generated matrices, eigenvectors, and perturbations reproducible.

The synthetic SPD matrices have the form

$$
A=Q\Lambda Q^T,
$$

where $Q$ is obtained from the QR factorization of a Gaussian random matrix and $\Lambda$ is prescribed.

Matrix-vector products are implemented using function handles, so the dense matrix $A$ does not need to be explicitly formed after $Q$ and the eigenvalues have been generated.

## Citation

If you use this code, please cite:

> H. Elzayyadi and J. M. Tabeart,  
> *The Influence of the Cluster Point on Rounding Errors and Sensitivity in the Spectral Limited-Memory Preconditioner.*

Publication details and DOI will be added when available.

## License

`pcg_modified.m` is derived from MATLAB's `pcg` implementation and retains the corresponding MathWorks copyright notice. The modification adds the output `x_hist`, which records the initial iterate and the CG iterate produced at each iteration and is used to compute relative forward-error histories.

The MATLAB helper routines used by `pcg_modified.m`, including `iterchk` and `iterapp`, are not distributed with this repository and are expected to be available through the user's MATLAB installation.

Users should consult the applicable MathWorks license terms regarding use and redistribution of MATLAB-derived code.
