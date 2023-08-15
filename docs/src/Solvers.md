# Solvers

Having defined a `SimpleLiquid` and a `Closure`, one should choose a method by which to solve the equations. The implemented methods are
`Exact`, `FourierIteration`, `NgIteration`. 

All solvers in some way, need to define a grid on which to solve the equations. This is done using the keyword arguments for `M` and `dr`, which represent respectively the number of grid points, and the spacing between them. Some solvers have additional settings, such as a tolerance. Using powers of 2 for `M` typically give the best performance (for 3D systems). For $N$ dimensional systems, the grid is constructed such that all `M` points lie below `M*dr`, but `dr` is not the exact grid spacing.

For the implemented cases, `Exact` solves the system exactly, or throws an error if the method is not implemented. 

```@docs
Exact
```

The methods `FourierIteration` and `NgIteration` both use recursive iteration to find improved estimates of the solution using Fourier Transforms. `NgIteration` uses a scheme to accelerate convergence, see Ref. [1]. 

```@docs
FourierIteration
NgIteration
```

## Meta-solvers

```@docs
DensityRamp
```