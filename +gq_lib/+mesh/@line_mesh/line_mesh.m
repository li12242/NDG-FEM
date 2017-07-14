classdef line_mesh < ndg_lib.mesh.line_mesh & gq_lib.mesh.mesh
    %LINE_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = line_mesh(cell, varargin)
            obj = obj@ndg_lib.mesh.line_mesh(cell, varargin{:});
            obj = obj@gq_lib.mesh.mesh(cell, varargin{:});
        end% func
    
    end
    
    methods(Hidden, Access=protected)
        function [nxq, nyq, nzq, wJsq] = face_scal(obj, vx, vy, vz, EToV)
            xb = obj.x(obj.cell.Fmask, :);
            nxq = ones(obj.cell.Nfptotal, obj.K);
            % Define outward normals
            [~, ind] = min(xb);
            nxq(ind, :) = -1;
            wJsq = ones(size(nxq));

            nyq = zeros(obj.cell.Nface, obj.K);
            nzq = zeros(obj.cell.Nface, obj.K);
        end% func
    end
    
end

