var documenterSearchIndex = {"docs":
[{"location":"HighDensities.html#Solving-at-high-densities","page":"Solving at high densities","title":"Solving at high densities","text":"","category":"section"},{"location":"ExtendingClosures.html#Defining-your-own-closure","page":"Defining your own closure","title":"Defining your own closure","text":"","category":"section"},{"location":"Potentials.html#Interaction-Potentials","page":"Interaction Potentials","title":"Interaction Potentials","text":"","category":"section"},{"location":"GeneralWorkflow.html#General-Workflow","page":"General Workflow","title":"General Workflow","text":"","category":"section"},{"location":"Solvers.html#Solvers","page":"Solvers","title":"Solvers","text":"","category":"section"},{"location":"ExtendingPotentials.html#Defining-your-own-potentials","page":"Defining your own potentials","title":"Defining your own potentials","text":"","category":"section"},{"location":"Closures.html#Closures","page":"Closures","title":"Closures","text":"","category":"section"},{"location":"SingleCompLJ.html#First-steps","page":"First steps","title":"First steps","text":"","category":"section"},{"location":"API.html#The-OrnsteinZernike-Module","page":"API","title":"The OrnsteinZernike Module","text":"","category":"section"},{"location":"API.html","page":"API","title":"API","text":"OrnsteinZernike","category":"page"},{"location":"API.html#OrnsteinZernike","page":"API","title":"OrnsteinZernike","text":"A generic solver package for Ornstein-Zernike equations from liquid state theory\n\n\n\n\n\n","category":"module"},{"location":"API.html#Module-Index","page":"API","title":"Module Index","text":"","category":"section"},{"location":"API.html","page":"API","title":"API","text":"Modules = [OrnsteinZernike]\nOrder   = [:constant, :type, :function, :macro]\nPrivate = false","category":"page"},{"location":"API.html#Detailed-API","page":"API","title":"Detailed API","text":"","category":"section"},{"location":"API.html","page":"API","title":"API","text":"Modules = [OrnsteinZernike]\nOrder   = [:constant, :type, :function, :macro]\nPrivate = false","category":"page"},{"location":"API.html#OrnsteinZernike.DensityRamp","page":"API","title":"OrnsteinZernike.DensityRamp","text":"DensityRamp <: Method\n\nSolves the system by iteratively solving systems of increasing density, using the previous solution as initial guess at a higher density. This may help deal with convergence issues\n\nArguments\n\nmethod: method by which to solve the system for individual densities.\ndensities: densities to consider. Must be a vector of increasing values.\nverbose: whether to print information.\n\nExample: DensityRamp(NgIteration(), [0.1, 0.3, 0.4]; verbose=false)\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.Exact","page":"API","title":"OrnsteinZernike.Exact","text":"Exact <: Method\n\nSolves the system exactly. This is only implemented for specific systems.\n\nConstruct using Exact(;  M=2^10, dr = sqrt(π/(M+1))/(2π)) Here, M is the number of points that the exact solution is evaluated on, and dr is the grid spacing. These are used to perform fourier transformations.\n\nExamples\n\nmethod = Exact() method = Exact(M=1000) method = Exact(M=1000, dr=0.01)\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.FourierIteration","page":"API","title":"OrnsteinZernike.FourierIteration","text":"FourierIteration <: Method\n\nSolves the system by recursive iteration in Fourier Space. Essentially, the algorithm is:\n\nguess an initial γ(r)\nfind c(r) using the closure relation\nfourier transform to get ĉ(k)\nfind γ(k) using the OZ-eq in k-space\ncompute γ(r) with a inverse fourier transform\ncompare with previous value, if not converged go to 2.\n\nOptionally, a mixing rule is used to mix the new and previous iteration of c(r) in step 2. \n\nArguments:\n\nM::Int: number of points discretize the solution on \ndr::Float64: grid spacing in real space\nmixing_parameter::Float64: mixing parameter for iteration mixing. A value of 1 is no mixing. Must be between 0 and 1. \nmax_iterations::Int64: maximal number of iterations \ntolerance::Float64: tolerance to be reached\nverbose::Bool: whether or not to print convergence information\n\ndefault: FourierIteration(; mixingparameter=0.5, maxiterations=10^5, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.HypernettedChain","page":"API","title":"OrnsteinZernike.HypernettedChain","text":"HypernettedChain\n\nImplements the Hypernetted Chain closure c(r) = (f(r)+1)*exp(γ(r)) - γ(r) - 1, or equivalently b(r) = 0.\n\nExample: closure = HypernettedChain()\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.MeanSphericalApproximation","page":"API","title":"OrnsteinZernike.MeanSphericalApproximation","text":"MeanSphericalApproximation\n\nImplements the MSA closure c(r) = -βu(r), or equivalently b(r) = ln(γ(r) - βu(r) + 1) - γ(r) +  βu(r).\n\nExample: closure = MeanSphericalApproximation()\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.MultiComponentHardSpheres","page":"API","title":"OrnsteinZernike.MultiComponentHardSpheres","text":"MultiComponentHardSpheres\n\nImplements the hard-sphere pair interaction uᵢⱼ(r) = inf for r < Dᵢⱼ and uᵢⱼ(r) = 0 for r > Dᵢⱼ for a multicomponent system.\n\nExpects a vector Dᵢ of diameters for each of the species. An additive mixing rule is used (Dᵢⱼ = (Dᵢ+Dⱼ)/2).\n\nExample: closure = MultiComponentHardSpheres([0.8, 0.9, 1.0])\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.NgIteration","page":"API","title":"OrnsteinZernike.NgIteration","text":"NgIteration <: Method\n\nSolves the system by recursive iteration in Fourier Space, and uses the Ng acceleration method. Essentially, the algorithm is:\n\nguess an initial γ(r)\nfind c(r) using the closure relation\nfourier transform to get ĉ(k)\nfind γ(k) using the OZ-eq in k-space\ncompute γ(r) with a inverse fourier transform\nuse Ng's method to provide a next guess for γ\ncompare with previous value, if not converged go to 2.\n\nArguments:\n\nM::Int: number of points discretize the solution on \ndr::Float64: grid spacing in real space\nN_stages::Int: Number of previous values to take into account for step 6. A higher number should lead to faster convergence, yet more computation time per iteration.\nmax_iterations::Int64: maximal number of iterations \ntolerance::Float64: tolerance to be reached\nverbose::Bool: whether or not to print convergence information\n\ndefault: NgIteration(; Nstages=3, maxiterations=10^3, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.OZSolution","page":"API","title":"OrnsteinZernike.OZSolution","text":"OZSolution\n\nHolds the solution of an Ornstein Zernike problem. \n\nFields:\n\nr: vector of distances\nk: vector of wave numbers\ngr: radial distribution function    \nSk: static structure factor\nck: direct correlation function in k space\ncr: direct correlation function in real space\n\nif the system was a single-component system, gr, Sk, ck and cr are vectors.  If instead the system was a multicomponent one, they are three dimensional vectors,  where the first dimension contains the values along r, and the second and third dimension contain the data for the species.\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.PercusYevick","page":"API","title":"OrnsteinZernike.PercusYevick","text":"PercusYevick\n\nImplements the Percus-Yevick closure c(r) = f(r)*(1+γ(r)), or equivalently b(r) = ln(1 + γ(r)) - γ(r).\n\nExample: closure = PercusYevick()\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.SimpleLiquid","page":"API","title":"OrnsteinZernike.SimpleLiquid","text":"SimpleLiquid{dims, ...} <: System\n\nHolds information about a homogeneous, isotropic system with radially symmetric pair interactions. dims is the dimensionality.\n\nConstruct using\n\nSimpleLiquid(dims, ρ, kBT, potential)\n\nFields:\n\nρ: number density, must be either a Number in case of a single component system, or a Vector in case of a mixture. In the latter case, each element contains the number density of the respective component.\nkBT: thermal energy\npotential::Potential: the interaction potential.  \n\nExamples:\n\nρ = 0.5; kBT = 1.1; dims = 3\npot = SingleComponentHardSpheres()\nsystem = SimpleLiquid(dims, ρ, kBT, pot)\n\nρ = [0.5, 0.1]; kBT = 5.2; dims = 3\npot = MultiComponentHardSpheres([1.0, 0.8])\nsystem = SimpleLiquid(dims, ρ, kBT, pot)\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.SingleComponentHardSpheres","page":"API","title":"OrnsteinZernike.SingleComponentHardSpheres","text":"SingleComponentHardSpheres\n\nImplements the hard-sphere pair interaction u(r) = inf for r < 1 and u(r) = 0 for r > 1.\n\nExample: closure = SingleComponentHardSpheres()\n\n\n\n\n\n","category":"type"},{"location":"API.html#OrnsteinZernike.compute_compressibility-Union{Tuple{P}, Tuple{T2}, Tuple{T1}, Tuple{dims}, Tuple{OZSolution, SimpleLiquid{dims, 1, T1, T2, P}}} where {dims, T1, T2, P}","page":"API","title":"OrnsteinZernike.compute_compressibility","text":"compute_compressibility(sol::OZSolution, system::SimpleLiquid)\n\nComputes the isothermal compressibility χ of the system\n\nuses the formula 1/χ = 1 - ρ ĉ(k=0) for single component systems and 1/χ = 1 - ρ Σᵢⱼ xᵢxⱼ ĉᵢⱼ(k=0) for mixtures. \n\n\n\n\n\n","category":"method"},{"location":"API.html#OrnsteinZernike.compute_excess_energy-Union{Tuple{P}, Tuple{T2}, Tuple{T1}, Tuple{OZSolution, SimpleLiquid{3, 1, T1, T2, P}}} where {T1, T2, P}","page":"API","title":"OrnsteinZernike.compute_excess_energy","text":"compute_excess_energy(sol::OZSolution, system::SimpleLiquid)\n\nComputes the excess energy per particle Eₓ, such that E = (dims/2kBT + Eₓ)N.\n\nuses the formula Eₓ = 1/2 ρ ∫dr g(r) u(r) for single component systems and Eₓ = 1/2 ρ Σᵢⱼ xᵢxⱼ ∫dr gᵢⱼ(r) uᵢⱼ(r) for mixtures. Here x is the concentration fraction ρᵢ/sum(ρ).\n\n\n\n\n\n","category":"method"},{"location":"API.html#OrnsteinZernike.compute_virial_pressure-Union{Tuple{P}, Tuple{T2}, Tuple{T1}, Tuple{OZSolution, SimpleLiquid{3, 1, T1, T2, P}}} where {T1, T2, P}","page":"API","title":"OrnsteinZernike.compute_virial_pressure","text":"compute_virial_pressure(sol::OZSolution, system::SimpleLiquid)\n\nComputes the pressure via the virial route\n\nuses the formula p = kBTρ - 1/6 ρ^2 ∫dr g(r) u'(r) for single component systems and p =  kBT Σᵢρᵢ - 1/6 Σᵢ ρᵢρⱼ ∫dr gᵢⱼ(r) u'ᵢⱼ(r) for mixtures.\n\n\n\n\n\n","category":"method"},{"location":"API.html#OrnsteinZernike.solve","page":"API","title":"OrnsteinZernike.solve","text":"solve(system::SimpleLiquid, closure::Closure, method::Method)\n\nSolves the system system using the closure closure with method method.\n\nsolve(system::SimpleLiquid, closure::Closure)\n\nSolves the system system using the closure closure with the default method NgIteration().\n\n\n\n\n\n","category":"function"},{"location":"HardSphereMixture.html#Hard-sphere-mixture","page":"Hard-sphere mixture","title":"Hard-sphere mixture","text":"","category":"section"},{"location":"index.html#OrnsteinZernike.jl","page":"Index","title":"OrnsteinZernike.jl","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"OrnsteinZernike.jl provides generic solvers for Ornstein-Zernike equations from liquid state theory.","category":"page"},{"location":"index.html#Installation","page":"Index","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"To install the package run ","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"import Pkg\nPkg.add(\"https://github.com/IlianPihlajamaa/OrnsteinZernike.jl\")","category":"page"}]
}
