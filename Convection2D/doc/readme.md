#Introduction

Two dimensional linear transport problem in order to assess the performance of the nodal DG method on triangular elements and on quadrilateral elements.

#Problem description

The governing equation

$$\frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u, x, t) = 0, \quad x\in\Omega$$

where *u* is the scalar variable, the flux term is

$$\mathbf{f} = \mathbf{a}(x) = (a_1(x)u, a_2(x)u).$$

1. advection velocity
$[a_1, a_2] = (-wy, wx)$, where $w = 5\pi/6$. The advection period is $T = \frac{2\pi}{w} = 2.4$.

2. initial condition
$u(\mathbf{x}, t = 0) = exp\left( -\sigma \left| \mathbf{x} - \mathbf{x}_c \right|^2 \right)$, where $\mathbf{x}_c = (0, 3/5)$, $\sigma = 125 \times 1000/33^2$.

3. computation domain
$[-L, L]$, where $L = 1$.

4. computation meshes
different resolutions: [20 x 20], [40 x 40], [60 x 60].
different mesh types:
	1. triangle mesh
	2. square mesh
	3. skewed-rectangular mesh

5. boundary condition
transmissive  boundary condition (zero gradient).


#Reference

[1] Wirasaet, D., Tanaka, S., Kubatko, E.J., Westerink, J.J., Dawson, C., 2010. A performance comparison of nodal discontinuous Galerkin methods on triangles and quadrilaterals. Int. J. Numer. Meth. Fluids 64, 1336–1362. doi:10.1002/fld.2376


