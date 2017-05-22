classdef detector1d < ndg_utility.detector.detector
    %DETECTOR1D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = detector1d(mesh, xd, dt, ftime, Nfield)
            yd = zeros(size(xd));
            zd = zeros(size(xd));
            obj = obj@ndg_utility.detector.detector(mesh, xd, yd, zd, dt, ftime, Nfield);
        end
        
        function [kd, rd, sd, td] = findlocate(obj)
            kd = ones(obj.Nd, 1);
            rd = ones(obj.Nd, 1);
            sd = ones(obj.Nd, 1);
            td = ones(obj.Nd, 1);
            for i = 1:obj.Nd
                dx = obj.mesh.x(obj.mesh.eidM) - obj.xd(i);
                kd(i) = find( dx(1,:).*dx(2,:) <= 0, 1 );
                vx = obj.mesh.vx( obj.mesh.EToV(:, kd(i)) );
                rd(i) = 2*( obj.xd(i) - 0.5*(vx(1) + vx(2)) )./( vx(2) - vx(1) );
            end
            
        end% func
        
    end
    
end

