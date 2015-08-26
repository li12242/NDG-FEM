classdef scalar_type < handle
    % Usage: s = scalar_type('scalar', mesh)
    properties
        val
        mesh
        name
        bc
    end %properties
    
    methods
        function obj = scalar_type(name, mesh)
            obj.name = name;
            obj.mesh = mesh;
            obj.bc = zeros(numel(mesh.mapB),1);
        end %function
        
        function setVal(obj, val)
            obj.val = val;
        end %function
        
        function faceVal = evaluateJump(obj)
            faceVal = obj.getInnerVal - obj.getOuterVal;
            faceVal(obj.mesh.mapB) = obj.getInnerVal(obj.mesh.mapB) - obj.bc(:);
        end
        
        function faceVal = evaluateMean(obj)
            faceVal = (obj.getInnerVal + obj.getOuterVal);
            faceVal(obj.mesh.mapB) = (obj.getInnerVal(obj.mesh.mapB) + obj.bc(:));
        end
        
        function obj = allocate(obj, size)
            obj.val = zeros(size);
        end %function
        
%         function obtainBC(obj, func)
%             obc.bc = feval(func);
%         end %function
    
        function val = getInnerVal(obj, ser)
            if nargin < 2
                val = obj.val(obj.mesh.vmapM);
            else 
                val = obj.val(obj.mesh.vmapM);
                val = val(ser);
            end
        end
        
        function val = getOuterVal(obj)
            val = obj.val(obj.mesh.vmapP);
        end
    end %methods
end %classdef