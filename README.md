# SciML
Workshop for SciML in Julia..

- Automated model discovery via universal differential equations.

# Universal approximation Theorem
Neural network are function expansions, fancy Taylor series
which are good for computing and bad for analysis.

> NN work well in high dimensions.

## Convolutional Neural Networks are structure assumptions.
- Universal Differential Equations
- Differential Equations defined in part by universal approximation
- Use all numerical features

## Toy model
> Paper: International Journal of Infections Diseases, Vol 93,
211-216

## SInDy - Sparse Identification of Dynamical Systems
> Paper: "Discovery governing equations from data by sparse
Identification of Nonlinear Dynamical Systems", PNAS 113.15
(2016): 3932-3937.

## ML Augmented Scientific Modeling
1. Identify known parts of a model, build UODE
2. Train a neural network to capture the missing mechanisms.
3. Sparse Identify the missing terms to mechanistic terms.
4. Verify the mechanisms are Scientific plausible
5. Extrapolate, do asymptotic analysis, predict bifurcations.
6. Get more data to verify the new terms.

### Examples
>Paper: Universal Battery Performance and Degradation Model for Electric Aircraft, DOI:10.26434/chemrxiv.12616169.v1
[https://chemrxiv.org/engage/chemrxiv/article-details/60c74d6cbdbb891026a3998f]

### Data driven quantification of Quarantine strength
>Paper: A machine learning aided global diagnostic and
comparative tool to assess effect of quarantine control in
COVID-19 spread.

## Beyond ODE
Universal Differential Algebraic Equations:
Encoding Physical Constraints

## Discretized PDE operators are convolutions

## Automatically Learning PDE from data: Universal PDE for Fisher-KPP

# Universal ODEs Accelerate Non-Newtonian Fluid Simulations

>Transform a systems of DAEs into parameterized system of ODEs 2x acceleration.

## Automated Climate Parameterizations
* Navier-Stokes are used in climate model.
* Attempt to solve this by "parameterizing" getting
a 1d approximation through averaging:
* Instead of picking w'c' replace it with NN and learn it
from small scale Simulations.

# UDE Methods Cover Accelerated Physics Informed Neural Networks Methods.
>Papers: "Multistep Neural Networks for Data-driven discovery of Nonlinear Dynamical Systems". Maziar Raissi,Paris Perdikaris,and George Em Karniadakis
"A comparative  study of physics-informed neural network models for learning unknown dynamics and constitutive relations" Ramakrishna Tipireddy, Paris Perdikaris, Pantos Stinis, and Alexandre Tartakovsky.

## Solving 1000 Dimensional Hamilton-Jacobi-Bellman via Universal SDEs

>Paper: Solving high-dimensional partial differential equations using deep learning [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6112690/]

UDE are a BLAS/LAPACK of SciML
Scientific Machine Learning requires efficient and accurate training of UDEs.

efficient and robust software for UDEs
# Packages
- DifferentialEquations.jl high Performance differential equations solvers
- DiffEqFlux.jl
- ModelingToolkit.jl
- NeuralPDE.jl: neural network solvers for PDEs, including automated physics-informed
- Catalyst.jl: high-performance differentiable modeling of chemical reaction networks.
- NBodySimulator.jl: high performance differential molecular dynamics
- DataDrivenDiffEq.jl: Koopman Dynamic mode decomposition (DMD)
and sparse identification (SInDy)

# SciML tools outperform ecosystems in high and low level languages

> Stiff2: [SciMLBenchmarks.jl]

Speed driver researchers to JuliaML.
Test problem Lorenz equations

* Stiff ODEs and DAEs
* Adjoint Methods
* Parallelism
* Event handling
* SDEs
* Delays

# SciML for flexibility, accuracy
## SciML tools do not rely on properties which can fail to hold

## Ill condition gradients cause difficulties in SciML

> Paper: Understanding and mitigating gradient pathologies
in physics informed neural networks. Sifan Wang, Yujun Teng, Paris Perdikaris.

> Note: Off the shelf ML tools will not work on stiff Scientific machine learning problems.

## DiffEqFlux has the features to handle stiff ill-conditioned Scientific problems.

- Rock Methods
- Implicit methods
- Multistep methods
- SSP methods (hyperbolic PDEs)
- Adaptive SDE solver (implicit higher order)
- Event handling

Neural ODE with batching on the GPU with high order adaptive
> Paper: Discretize-Optimize vs. Optimize-Discretize for Time series and Continuous Normalizing Flows. Derek Onken, Lars Ruthotto.
