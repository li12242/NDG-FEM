classdef Visual2d < AbstractVisual
    properties (SetAccess = protected)
        %> vertex list in sub-element
        SEToV
        %> num of sub-cell
        Ncell
        %> num of tri face
        Ntri
        %> m-by-3 face matrix
        tri
    end

    methods (Access = public)
        function obj = Visual2d(output)
            obj = obj@AbstractVisual(output);
            obj = obj.InitVisual2d(  );
        end

        drawResult( obj, step, fieldId );
        drawMesh( obj );
    end

    methods (Access = protected)
        % init
        obj = InitVisual2d( obj );
    end
end