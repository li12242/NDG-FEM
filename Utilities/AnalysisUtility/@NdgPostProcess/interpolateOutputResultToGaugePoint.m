%> \brief Return interpolation of numerical result at input gauge points.
%> \details Function `interpolateOutputResultToGaugePoint` interpolate the
%> numerical results at gauge points with three steps:
%> 1. The localtion of cell index and local coordinate [r, s] are
%> deternmined;
%> 2. The interpolation matrix Vg is calculated with its entries given by
%> \f$ [V_{g}]_{i,j} = \varphi_j(r_i, s_i), \; i = 1, \cdots, N_g \f$,
%> where \f$ N_g \f$ is the total number of gauge points.
%> 3. The interpolated results at gauge point \f$(x_i, y_i)\f$ are obtained 
%> with \f$ \sum_{j = 1}^{N_p} [V_{g}]_{i,j} \mathbf{U}(j, k_i) \f$, where
%> \f$ k_i \f$ is the cell index.
function [ gaugeValue ] = interpolateOutputResultToGaugePoint( obj, xg, yg, zg )
Ntime = obj.Nt;
Ng = numel( xg );
gaugeValue = zeros(Ntime, obj.Nvar, Ng);
for m = 1:obj.Nmesh
    [ cellId, Vg ] = obj.meshUnion.accessGaugePointLocation( xg, yg, zg );
    % get the output file name
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    % loop over all the output step
    for t = 1:Ntime
        fresult = ncread(obj.outputFile{m}, 'fphys', ...
            [1, 1, 1, t], [Np, K, obj.Nvar, 1]);
        for n = 1:Ng
            if (cellId(n) == 0)
                % the gauge point locates out of the mesh domain
                continue;
            end
            for fld = 1:obj.Nvar
                gaugeValue(t, fld, n) = Vg(n, :) * fresult(:, cellId(n), fld);
            end
        end
    end
end
end% func