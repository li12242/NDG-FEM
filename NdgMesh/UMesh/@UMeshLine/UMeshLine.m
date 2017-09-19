%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshLine < UMeshUnion

    methods(Hidden, Access = protected)
        function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = getElementalNodeInfo(obj)
            Np = obj.cell.Np; K = obj.K;
            xr = obj.cell.Dr*obj.x;
            J = xr; 
            rx = 1./xr;
            ry = zeros(Np, K);
            rz = zeros(Np, K);
            sx = zeros(Np, K);
            sy = zeros(Np, K);
            sz = zeros(Np, K);
            tx = zeros(Np, K);
            ty = zeros(Np, K);
            tz = zeros(Np, K);
        end
%         [nx, ny, nz, Js] = getElementalSurfaceInfo(obj, vx, vy, vz, EToV)
        
%         function [ EToFG ] = getElementalAdjacentFaceGlobalIndex(obj)
%             EToFG = zeros(obj.cell.Nface, obj.K);
%             for f = 1:obj.cell.Nface
%                 EToFG(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
%             end
%         end% func
    end% methods
    
    methods
        obj = refine(obj, refine_level);
        
        function obj = UMeshLine(cell, varargin)
            [cell, Nv, vx, K, EToV, EToR, EToBS] = readInput(cell, varargin{:});
            
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros
            
            obj = obj@UMeshUnion(cell, Nv, vx, vy, vz, K, EToV, EToR, EToBS);
        end% func
        
        function draw(obj)
            plot(obj.x, zeros(obj.cell.Np, obj.K), '.-');
        end
    end
    
end

function [cell, Nv, vx, K, EToV, EToR, EToBS] = readInput(cell, varargin)
msgID = 'UMeshLine:Inputerror';
% check cell type
if (~isa(cell, 'StdLine'))
    msgtext = 'The input cell should be StdLine object.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if(nargin == 2) % case input
    casename = varargin{1};
    [Nv, vx, K, EToV, EToR, EToBS] = read_from_file(casename);
elseif(nargin == 7) % parameters input
    [Nv, vx, K, EToV, EToR, EToBS] = checkInputVars(varargin{:});
else
    msgtext = 'The number of input should be 2 or 7.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

function [Nv, vx, K, EToV, EToR, EToBS] = checkInputVars(varargin)
msgID = 'UMeshLine:Inputerror';

Nv = varargin{1};
vx = varargin{2};
if (numel(vx) ~= Nv)
    msgtext = 'The number of vertices in vx is not equal to Nv.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

K = varargin{3};
EToV = varargin{4};
EToR = varargin{5};
EToBS = varargin{6};

if ( size(EToV, 2) ~= K )
    msgtext = 'The number of elements in EToV is not equal to K.';
    ME = MException(msgID, msgtext);
    throw(ME);
elseif (numel(EToR, 2) ~= K)
    msgtext = 'The number of elements in EToR is not equal to K.';
    ME = MException(msgID, msgtext);
    throw(ME);
elseif (numel(EToBS, 2) ~= K) 
    msgtext = 'The number of elements in EToBS is not equal to K.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func


function [Nv, vx, K, EToV, EToR, EToBS] = read_from_file(casename)
end
