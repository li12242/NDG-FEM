classdef NdgConDetector < NdgAbstractDetector
    
    methods
        function flg = assembleTroubleCell( obj, fphys, fldId )  
            flg = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                flag{m} = zeros();
                
                flg{m} = uint( flag{m} );
            end
            
        end
    end
    
end

