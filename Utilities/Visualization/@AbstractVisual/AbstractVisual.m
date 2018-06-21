classdef AbstractVisual < handle
    properties (SetAccess = protected)
        mesh
        outputObj
        drawHandle
    end

    methods (Abstract, Access = public)
        drawResult( obj, field );
        drawOutputResult( obj, outputObj, step, fieldId );
        drawMesh( obj );
    end% methods

    methods (Access = public)
        function obj = AbstractVisual( mesh )
            obj.mesh = mesh;
        end
        
        function drawOutputResultAll( obj, outputObj, stepInterval, fieldId )
            Nstep = outputObj.outputStep;
            for n = 1 : stepInterval : Nstep
                obj.drawOutputResult( outputObj, n, fieldId );
                drawnow;
            end
        end
    end% methods
end