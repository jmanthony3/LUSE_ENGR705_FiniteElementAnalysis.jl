# LUSE_ENGR705_FiniteElementAnalysis.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmanthony3.github.io/LUSE_ENGR705_FiniteElementAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmanthony3.github.io/LUSE_ENGR705_FiniteElementAnalysis.jl/dev/)
[![Build Status](https://github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/jmanthony3/LUSE_ENGR705_FiniteElementAnalysis.jl/actions/workflows/Documentation.yml)

This package provides:
- Command line interface for simple, two-dimensional, solid mechanics Finite Element Analysis (FEA) problems.
- Several functions with method dispatching according to FEA element type.
- Support for "line" and "beam" element types.
- Customization for which LinearSolve.jl [algorithm](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/) to employ.

Key comments:
- Currently, the global stiffness matrix, $\mathbf{K}$, is constructed according to the *Minimum Potential Energy Method*.
- The `solve()` function heavily leverages the [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) package to solve the linear system of equations: $\vec{F} = \mathbf{K}\vec{U}$.
- Package can currently solve most truss analysis problems and some beams with distributed loading.

## Roadmap
- [x] 2D
  - [x] "line" element types
  - [x] "beam" element types
  - [ ] "triangular" element types
  - [ ] Checks for symmetry
- [ ] 3D
  - [ ] Extend 2D elements for third component direction
  - [ ] Other 3D element types:
    - [ ] "quadralateral"
    - [ ] "tetrahedronal"
    - [ ] etcetera...
- [ ] Condense `prepareelements_*` functions to single function with method dispatching for keyword assignment `elementtypes="line"`
- [ ] Shear and moment diagrams
- [ ] Support for mixed element types in `solve()`
- [ ] Mesh refinement function
- [ ] Option to solve FEA problems according to other techniques
- [ ] Extensive exception handling and error testing for each method of package

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
