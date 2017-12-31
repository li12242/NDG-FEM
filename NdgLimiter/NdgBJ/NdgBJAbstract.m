classdef NdgBJAbstract < NdgAbstractLimiter
    
    properties
        K
    end
    
    methods
        function obj = NdgBJAbstract( mesh )
            obj = obj@NdgAbstractLimiter( mesh );
        end
    end
    
    methods( Hidden )
        function [ cmean, cmax, cmin ] = accessCellExtrema( obj, fphys, fldId )
            cmax = cell( obj.Nmesh, 1 );
            cmin = cell( obj.Nmesh, 1 );
            
            cmean = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                cmean{m} = obj.meshUnion(m).GetMeshAverageValue...
                    ( fphys{m}(:,:,fldId) );
            end
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                hca = zeros(mesh.cell.Nface, mesh.K);
                ind = ( mesh.EToM == m );
                hca( ind ) = cmean{m}( mesh.EToE( ind ) );
                for e = 1:numel(mesh.edgeUnion)
                    edge = mesh.edgeUnion(e);
                    m2 = edge.FToM;
                    for n = 1:edge.M
                        k1 = edge.FToE(1, n);
                        k2 = edge.FToE(2, n);
                        f1 = edge.FToF(1, n);
                        f2 = edge.FToF(1, n);
                        
                        hca(f1, k1) = cmean{m2}(f2, k2);
                    end
                end
                
                cmax{m} = max( [cmean{m}; hca] );
                cmin{m} = min( [cmean{m}; hca] );
            end
        end
    end
    
end
