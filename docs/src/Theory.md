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

$h_{ij}(r) = c_{ij}(r) + \sum_l \rho_l \int d\textbf{r}' c_{il}(\textbf{r}')h_{lj}(|\textbf{r} - \textbf{r}'|) $

links the pair correlation function $g_{ij}(r) = h_{ij}(r)+1$ to the species specific number density $\rho_{i}$ and the direct correlation function $c_{ij}(r)$. In order to solve it, a second relation between $h_{ij}(r)$ and $c_{ij}(r)$ must be specified. This is called the closure relation. 

In practise, the closure relation typically takes the form $c_{ij}(r) = f(\gamma_{ij}(r), r, u_{ij}(r))$, where $u_{ij}(r)$ is the interaciton potential and $\gamma_{ij}(r) = h_{ij}(r) - c_{ij}(r)$ is the indirect correlation function. Sometimes, instead the closure is defined for the bridge function $b(r)$, and the closure relation is then given by $h_{ij}(r) - 1 = \exp\left(-\beta u_{ij}(r) + h_{ij}(r) - c_{ij}(r) + b_{ij}(r) \right)$, where $\beta = 1/k_BT$. 

## Fourier Transforms

## Thermodynamic properties

