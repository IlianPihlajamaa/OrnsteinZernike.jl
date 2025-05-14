# Solvers

Having defined a `SimpleLiquid` and a `Closure`, one may choose a method by which to solve the equations. The implemented methods are
`Exact`, `FourierIteration`, `NgIteration`. If no method is given the `solve` function, it will use the default `NgIteration`.

All solvers in some way need to define a grid on which to solve the equations. This is done using the keyword arguments for `M` and `dr`, which represent respectively the number of grid points, and the spacing between them. Some solvers have additional settings, such as a tolerance. It is important to ensure that `dr` is small compared to the features of the interaction potential, and that `M*dr` is sufficiently large. For $N$ dimensional systems, the grid is constructed such that all `M` points lie below `M*dr`, but `dr` may not be the exact grid spacing between all points.

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
TemperatureRamp
```

## References
[1] Ng, K. C. (1974). Hypernetted chain solutions for the classical one‐component plasma up to Γ= 7000. The Journal of Chemical Physics, 61(7), 2680-2689.