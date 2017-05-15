# Wet/Dry treatment

## Identify the refinement element

The following process is in function

	refineflag = TransiteCellIdentify(mesh, h, bedElva)

where the variable `refineflag` will return a bool vector flag for wet/dry refinement element.

###1. identify the transitation element
For define the transitation (wet/dry interface) element, the element must satisfy the following condition

$$\left\{ \begin{matrix}
h_i > h_{dry} \quad \text{for} \quad \forall x_i \in \Omega_i \cr
\Omega_{i-1} \in \text{wet cell} \quad \& \quad \Omega_{i+1} \in \text{dry cell}
\end{matrix}\right.
$$

or

$$
\left\{ \begin{matrix}
h_i > h_{dry} \quad \text{for} \quad \forall x_i \in \Omega_i \cr
\Omega_{i-1} \in \text{dry cell} \quad \& \quad \Omega_{i+1} \in \text{wet cell}
\end{matrix}\right.
$$

The first condition requires that there is water in the transitation element $\Omega_i$. The second condition require the adjacent element has different wet-dry status.

### 2. the refinement condition

For the transitation element which need refinement, the water depth satisfy the following condition

$$ \bar{h}_i + B_{i+1/2} < \mathrm{max}(B_i, B_{i+1}) $$

## h refinement

Function

	[new_mesh, h, q, localEleIndex] = Hrefine1D(refineflag, mesh, h, q)

refine the transitation and returns the new mesh object with variable `h` & `q`. The variable `localEleIndex` is a matrix store the information of local refinement.

### 1. get new vertex index, location and EToV

### 2. get new EToE and EToF

### 3. interpolate h and q to new local refinement element


## Combine the local refined element
