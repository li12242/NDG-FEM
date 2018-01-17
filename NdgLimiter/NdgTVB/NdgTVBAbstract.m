classdef NdgTVBAbstract < NdgAbstractLimiter
    %NDGTVBABSTRACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M
    end
    
    methods
        function obj = NdgTVBAbstract( mesh, M )
            obj = obj@NdgAbstractLimiter(mesh);
            obj.M = M;
        end
    end
    
    methods( Access = protected )
        function [ hc, hca ] = matEvaluateIntegralMean( obj, fphys, fldId )
            hc = cell( obj.Nmesh );
            hca = cell( obj.Nmesh );
            for m = 1:obj.Nmesh
                hc{m} = obj.meshUnion(m).GetMeshAverageValue...
                    ( fphys{m}(:,:,fldId) );
            end
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                hca{m} = zeros(mesh.cell.Nface, mesh.K);
                ind = ( mesh.EToM == m );
                hca{m}( ind ) = hc{m}( mesh.EToE( ind ) );
                for e = 1:numel(mesh.edgeUnion)
                    edge = mesh.edgeUnion(e);
                    m2 = edge.FToM;
                    for n = 1:edge.M
                        k1 = edge.FToE(1, n);
                        k2 = edge.FToE(2, n);
                        f1 = edge.FToF(1, n);
                        %f2 = edge.FToF(1, n);
                        
                        hca{m}(f1, k1) = hc{m2}(k2);
                    end
                end
            end
        end% func
    end
    
    methods( Access = protected, Static )
        function [ m ] = minmod( v )
            num = size(v,1);
            m = zeros(1, size(v,2));
            s = sum(sign(v), 1)/num ;
            
            ids = ( abs(abs(s)-1)<1e-10 );
            tmp = s.*min( abs(v), [], 1);
            m(ids) = tmp(ids);
        end
    end
    
end

