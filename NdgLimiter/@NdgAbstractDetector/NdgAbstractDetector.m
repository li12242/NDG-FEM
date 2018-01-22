classdef NdgAbstractDetector < handle
    
    properties
        meshUnion
        Nmesh
    end
    
    methods
        function obj = NdgAbstractDetector( meshUnion )
            obj.Nmesh = numel( meshUnion );
            obj.meshUnion = meshUnion;
        end
    end
    
    methods( Abstract )
        flg = assembleTroubleCell( obj, fphys, fldId )            
    end
    
end

