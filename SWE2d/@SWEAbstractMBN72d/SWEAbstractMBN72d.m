classdef SWEAbstractMBN72d < SWEAbstractDB2d
    
    properties(Constant)
        %> Physical field - {h, hu, hv, b, bx, by, eta}
        Nfield = 7
    end
    
    methods
        function obj = SWEAbstractMBN72d( )
            obj = obj@SWEAbstractDB2d( );
        end
    end
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            ind = (mesh.EToR == int8( NdgRegionType.Wet ));
            E = zeros( mesh.cell.Np, mesh.K, obj.Nvar );
            G = zeros( mesh.cell.Np, mesh.K, obj.Nvar );

            h2 = 0.5 * obj.gra * ( fphys(:, ind, 1).^2 - fphys(:, ind, 4).^2 );
            huv = fphys(:, ind, 2) .* fphys(:, ind, 3) ./ fphys(:, ind, 1);
            E(:, ind, 1) = fphys(:, ind, 2);
            G(:, ind, 1) = fphys(:, ind, 3);
            E(:, ind, 2) = fphys(:, ind, 2).^2./fphys(:, ind, 1) + h2;
            G(:, ind, 2) = huv;
            E(:, ind, 3) = huv;
            G(:, ind, 3) = fphys(:, ind, 3).^2./fphys(:, ind, 1) + h2;
        end

%         [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fphys, fext );
%         [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fphys );
    end
    
    methods(Access=protected)
        function fphys = matEvaluateLimiter( obj, fphys )
            fphys = obj.limiter.matLimit( fphys, 1 );
            fphys = obj.limiter.matLimit( fphys, 2 );
            fphys = obj.limiter.matLimit( fphys, 3 );
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,7) = fphys{m}(:,:,1) + fphys{m}(:,:,4);
            end
            fphys = obj.limiter.matLimit( fphys, 7 ); % enforce the elevation

            for m = 1:obj.Nmesh % set the new bottom topography
                mesh = obj.meshUnion(m);
                ind = (mesh.EToR == int8( NdgRegionType.Wet ));
                fphys{m}(:,ind,4) = fphys{m}(:,ind,7) - fphys{m}(:,ind,1);
                
                fphys{m}(:,ind,5) = ...
                    mesh.rx(:, ind) .* (mesh.cell.Dr * obj.fphys{m}(:,ind,4)) + ...
                    mesh.sx(:, ind) .* (mesh.cell.Ds * obj.fphys{m}(:,ind,4));
                fphys{m}(:,ind,6) = ...
                    mesh.ry(:, ind) .* (mesh.cell.Dr * obj.fphys{m}(:,ind,4)) + ...
                    mesh.sy(:, ind) .* (mesh.cell.Ds * obj.fphys{m}(:,ind,4));
            end
        end
    end
    
end

