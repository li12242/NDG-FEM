classdef AbstractVisual < handle
    properties (SetAccess = protected)
        mesh
        outputObj
        drawHandle
    end

    methods (Abstract, Access = public)
        drawResult( obj, field );
        drawOutputResult( obj, step, fieldId );
        drawMesh( obj );
    end% methods

    methods (Access = public)
        function obj = AbstractVisual( output )
            obj.outputObj = output;
            obj.mesh = output.mesh;
        end
        
        function drawOutputResultAll( obj, stepInterval, fieldId )
            Nstep = obj.outputObj.outputStep;
            for n = 1 : stepInterval : Nstep
                obj.drawOutputResult( n, fieldId );
                drawnow;
            end
        end
    end% methods
end