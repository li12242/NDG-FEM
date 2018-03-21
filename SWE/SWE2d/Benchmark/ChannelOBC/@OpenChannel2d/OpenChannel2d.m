classdef OpenChannel2d < SWEWD2d
    
    properties(Constant)
        hmin = 1e-4;
        gra = 9.81;
        H = 40; % mean depth
        eta = 0.2; % amplitude
        L = 2e3;
        W = 500;
        % mesh center and rotation angle
        theta = 0;
        xc = 0;
        yc = 0;
    end
    
    methods
        function obj = OpenChannel2d( N, type )
            obj = obj@SWEWD2d();
            mesh = makeUniformMesh( N, type,  )
        end
    end
    
end

