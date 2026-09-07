# Spectral-LMP Numerical Experiments

MATLAB code accompanying numerical experiments on the influence of the cluster
point in the **Spectral Limited-Memory Preconditioner (Spectral-LMP)**.

The experiments investigate:

- finite-precision effects in the application of Spectral-LMP;
- sensitivity to inaccurate dominant eigenvectors;
- the effect of perturbation magnitude and direction;
- several choices of the cluster point $\theta$;
- approximate spectral information computed using randomized Nyström and
  REVD-ritzit methods.

## Repository structure

### Core Spectral-LMP routines

`Spectral_LMP.m`  
Applies the Spectral-LMP preconditioner using the supplied dominant
eigenvectors, eigenvalues, and cluster point.

`run_LMP_experiment.m`  
Runs CG/PCG for several choices of the cluster point, including

- $\theta_{WM,e}$;
- $\theta_{WM,1/\lambda}$;
- $\lambda_k$;
- $\lambda_{k+1}$;
- $\theta_r$;
- $\theta_m$.

`pcg_modified.m`  
Modified version of MATLAB's built-in `pcg` routine. It is retained here so
that the numerical experiments use the same MATLAB PCG implementation employed
in the study. The modification adds the output `x_hist`, which stores the
initial iterate and the CG iterate produced at each iteration. This history is
used to compute the relative forward error throughout the iteration.

The routine relies on helper functions used by MATLAB's original `pcg`
implementation, including routines such as `iterchk`, `iterapp`, and
`itermsg`. These helper routines are not included in this repository and are
expected to be available through the user's MATLAB installation.

This file is derived from MathWorks code and retains the corresponding
copyright notice. It is included to document and reproduce the computational
setup used for the numerical experiments. Users should consult the applicable
MathWorks license terms regarding use and redistribution.

### Plotting routines

`plot_convergence.m`  
Plots the relative forward-error histories for the different cluster-point
choices.

`plot_eigenvector_diagnostics.m`  
Visualizes the accuracy and subspace structure of the approximate dominant
eigenvectors.

The routine produces heat maps of

$$
|\widehat V_k^T Q_k|
$$

and

$$
|\widehat V_k^T Q_{k+1:n}|,
$$

which illustrate, respectively, mixing within the dominant eigenspace and
leakage into the complement eigenspace.

It also plots, for each approximate eigenvector $\widehat v_i$, the total
perturbation magnitude

$$
\tau_i = \lVert \widehat v_i-v_i\rVert _2,
$$

together with its dominant- and complement-subspace components,

$$
\eta_i = \lVert Q_k^T(\widehat v_i-v_i)\rVert_2,
$$

and

$$
\delta_i = \lVert Q_{k+1:n}^T(\widehat v_i-v_i)\rVert_2.
$$

These diagnostics distinguish the overall perturbation magnitude from its
direction relative to the dominant and complement eigenspaces.

## Experiments

### `ExactLMPvsLMP.m`

Compares

1. unpreconditioned CG;
2. CG applied to an exactly clustered system; and
3. PCG using the practical Spectral-LMP.

The exactly clustered system and the Spectral-LMP-preconditioned system have
the same clustered spectrum in exact arithmetic. Their different numerical
behavior therefore illustrates rounding-error effects in the practical
application of the preconditioner.

Run with

```matlab
ExactLMPvsLMP
```

### `LMPRoundingError.m`

Investigates rounding errors in the dominant and complement subspaces and
compares observed errors with the corresponding theoretical bounds.

Run with

```matlab
LMPRoundingError
```

### `SyncSLMP_CG.m`

Runs synthetic experiments with perturbed dominant eigenvectors.

The script

1. generates random perturbation directions;
2. controls the complement-subspace component of the perturbations;
3. scales the perturbations according to prescribed values $\tau_i$;
4. reorthogonalizes the perturbed eigenvectors;
5. computes the actual perturbation magnitudes after reorthogonalization;
6. runs the Spectral-LMP/CG experiments; and
7. computes dominant-subspace mixing and complement-subspace leakage
   diagnostics.

Run with

```matlab
SyncSLMP_CG
```

The script runs all three synthetic perturbation experiments defined in
`SyncSLMP_CG.m` and produces the corresponding convergence and eigenvector
diagnostics.

### `Nystrom_Experiment.m`

Computes approximate dominant eigenpairs using the randomized Nyström method
and investigates the corresponding Spectral-LMP convergence and eigenvector
diagnostics.

Run with

```matlab
Nystrom_Experiment
```

The randomized eigensolver is implemented in

```text
randomized_nystrom.m
```

### `Ritzit_Experiment.m`

Computes approximate dominant eigenpairs using the REVD-ritzit method and
performs the same convergence and eigenspace diagnostics.

Run with

```matlab
Ritzit_Experiment
```

The randomized eigensolver is implemented in

```text
REVD_ritzit.m
```

## Requirements

The code is written in MATLAB.

The experiments use standard dense linear-algebra routines including

```text
qr
svd
eig
chol
```

`LMPRoundingError.m` currently uses `normrnd`. This function belongs to the
Statistics and Machine Learning Toolbox. If desired, the relevant normally
distributed random vectors can instead be generated using `randn`.

The supplied experiments use matrices of dimension $n=5000$. Since a dense
$5000\times5000$ orthogonal matrix is generated, the experiments can require
substantial memory and runtime.

## Running the code

Place all `.m` files in the same directory and make that directory the MATLAB
working directory, or add it to the MATLAB path.

For example,

```matlab
cd path/to/Spectral-LMP
addpath(pwd)
```

Then run one of the experiment scripts, for example

```matlab
SyncSLMP_CG
```

or

```matlab
Nystrom_Experiment
```

## Reproducibility

The experiment scripts use fixed random-number-generator seeds of the form

```matlab
rng(...,'twister')
```

so that the randomly generated matrices, eigenvectors, and perturbations are
reproducible.

The synthetic SPD matrices have the form

$$
A = Q\Lambda Q^T,
$$

where $Q$ is obtained from the QR factorization of a Gaussian random matrix
and $\Lambda$ is prescribed. Matrix-vector products are performed using
function handles, so the dense matrix $A$ does not need to be formed
explicitly after $Q$ and the eigenvalues have been generated.

## Eigenvector perturbation diagnostics

The experiments distinguish between the magnitude and direction of the
eigenvector perturbations. For each approximate eigenvector, $\tau_i$ measures
the total perturbation magnitude, while $\eta_i$ and $\delta_i$ measure its
components in the dominant and complementary eigenspaces, respectively.

The heat maps of $|\widehat V_k^T Q_k|$ and
$|\widehat V_k^T Q_{k+1:n}|$ provide additional information about mixing
within the dominant eigenspace and leakage into the complementary eigenspace.
These diagnostics are used to relate the structure of the approximate
spectral information to the observed Spectral-LMP convergence.

## Citation

If you use this code, please cite the accompanying paper:

> H. Elzayyadi and J. M. Tabeart,  
> *The Influence of the Cluster Point on Rounding Errors and Sensitivity in the
> Spectral Limited-Memory Preconditioner.*

Publication details and DOI can be added here when available.

## License

A license for the original code in this repository should be added before
public release.

`pcg_modified.m` is derived from MATLAB's `pcg` implementation and is therefore
not presented as original code of the authors. It retains the MathWorks
copyright notice and is included because the experiments were performed using
this modified MATLAB implementation. The modification records the history of
the CG iterates through the additional `x_hist` output.

The MATLAB helper routines called by `pcg_modified.m` are not distributed with
this repository. Users should rely on their MATLAB installation for those
dependencies and should consult the applicable MathWorks license terms for the
MATLAB-derived code.
