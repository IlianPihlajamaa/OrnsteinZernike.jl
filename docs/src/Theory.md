# Theory

In this section, we describe our conventions and notation.

## Single component systems

The Ornstein zernike equation 
$h(r) = c(r) + \rho \int d\textbf{r}' c(\textbf{r}')h(|\textbf{r} - \textbf{r}'|) $
links the pair correlation function $g(r) = h(r)+1$ to the number density $\rho$ and the direct correlation function $c(r)$. In order to solve it, a second relation between $h(r)$ and $c(r)$ must be specified. This is called the closure relation. 

In practise, the closure relation typically takes the form $c(r) = f(\gamma(r), r, u(r))$, where $u(r)$ is the interaciton potential and $\gamma(r) = h(r) - c(r)$ is the indirect correlation function. Sometimes, instead the closure is defined for the bridge function $b(r)$, and the closure relation is then given by $h(r) - 1 = \exp\left(-\beta u(r) + h(r) - c(r) + b(r) \right)$, where $\beta = 1/k_BT$. 

## Mixtures
Everything above generalizes to the mixture case:

The Ornstein zernike equation 

$h_{ij}(r) = c_{ij}(r) + \sum_l \rho_l \int d\textbf{r}' c_{il}(\textbf{r}')h_{lj}(|\textbf{r} - \textbf{r}'|)$

links the pair correlation function $g_{ij}(r) = h_{ij}(r)+1$ to the species specific number density $\rho_{i}$ and the direct correlation function $c_{ij}(r)$. In order to solve it, a second relation between $h_{ij}(r)$ and $c_{ij}(r)$ must be specified. This is called the closure relation. 

In practise, the closure relation typically takes the form $c_{ij}(r) = f(\gamma_{ij}(r), r, u_{ij}(r))$, where $u_{ij}(r)$ is the interaciton potential and $\gamma_{ij}(r) = h_{ij}(r) - c_{ij}(r)$ is the indirect correlation function. Sometimes, instead the closure is defined for the bridge function $b(r)$, and the closure relation is then given by $h_{ij}(r) - 1 = \exp\left(-\beta u_{ij}(r) + h_{ij}(r) - c_{ij}(r) + b_{ij}(r) \right)$, where $\beta = 1/k_BT$. 

## Fourier Transforms

The ability to numerically solve the Ornstein-Zernike equation relies heavily on doing repeated Fourier Transforms. In arbitrary dimensions, these Fourier transforms can be written as Hankel transforms in the case that the argument is a radially symmetric function. In particular, in $d$ dimensions, the radial Fourier transform and its inverse are given by

$\hat{F}(k) = (2\pi)^{d/2} k ^{1-d/2}\int_0^\infty dr r^{d/2}J_{d/2-1}(kr)F(r)$

$F(r) = (2\pi)^{-d/2} r ^{1-d/2}\int_0^\infty dk k^{d/2}J_{d/2-1}(kr)\hat{F}(k),$

in which $J_{d/2-1}(x)$ is the bessel function of order $d/2-1$. In the special cases of 1 and 3 dimensions, the transform simplifies into:

$\hat{F}(k) = 2\int_0^\infty dr \cos(kr)F(r)$

$F(r) = \frac{1}{\pi}\int_0^\infty dk \cos(kr)\hat{F}(k),$

in 1$d$, and 

$\hat{F}(k) = \frac{4\pi}{k}\int_0^\infty dr r \sin(kr)F(r)$

$F(r) = \frac{1}{2\pi^2r}\int_0^\infty dk k\sin(kr)\hat{F}(k),$

in 3$d$. This package uses discrete versions of all of the above. 

## Thermodynamic properties

