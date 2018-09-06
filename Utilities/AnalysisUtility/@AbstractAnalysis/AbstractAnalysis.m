classdef AbstractAnalysis < handle
    %ABSTRACTPOSTPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        phys
    end
    
    properties ( SetAccess = protected )
        %> num of gauge points
        Ng
        %> num of mesh
        Nmesh
        %> num of physical field
        Nfield
        %> gauge point position
        xg, yg, zg
        %> mesh id of gauge point location
        gaugeMesh
        %> cell id of gauge point location
        gaugeCell
        %> interp matrix
        Vg
        %> number of gauge points in each mesh
        NgaugePerMesh
        %> gauge points index in each mesh
        gaugeIdPerMesh
    end

    methods (Abstract, Access = protected)
        [ gaugeMesh, gaugeCell ] = accessGaugePointLocation( obj );
        [ Vg ] = assessGaugeInterpMatrix( obj );
    end
    
    methods(Access = protected)
        function assessGaugePerMesh( obj )
            obj.NgaugePerMesh = zeros( obj.Nmesh, 1 );
            obj.gaugeIdPerMesh = zeros( obj.Ng, obj.Nmesh );
            for n = 1:obj.Ng
                meshId = obj.gaugeMesh(n);
                obj.NgaugePerMesh( meshId ) = obj.NgaugePerMesh( meshId ) + 1;
                obj.gaugeIdPerMesh( obj.NgaugePerMesh( meshId ) ) = n;
            end            
        end
    end
    
    methods
        function obj = AbstractAnalysis( phys )
            obj.phys = phys;
            obj.Nfield = phys.Nfield;
            obj.Nmesh = phys.Nmesh;
        end

        data = InterpGaugeResult( obj, fphys )
        data = InterpOutputResult2GaugePoint( obj, step )
    end
    
end

