classdef AdvRotationShockUniformMesh2d < AdvRotationUniformMesh2d
    
    methods
        function obj = AdvRotationShockUniformMesh2d(N, M, cellType)
            obj = obj@AdvRotationUniformMesh2d(N, M, cellType);
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fphys{m} = getExtFunc(obj, obj.meshUnion(m), 0);
            end
        end
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 2.4;
            option('LimiterType') = NdgLimiterType.Vert;
            option('temporalDiscreteType') = NdgTemporalIntervalType.Constant;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = 2/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('timeInterval') ...
                = sqrt(2)/obj.M/obj.w/(2*obj.N + 1);
            option('equationType') = NdgDiscreteEquationType.Weak;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
        
        
        function [ x, y, r ] = getRelativeDistanceToExactPosition( obj, mesh, theta, time )
            theta = theta + obj.w*time;
            x = obj.x0 + obj.rd*cos(theta);
            y = obj.y0 + obj.rd*sin(theta);
            r  = sqrt( (mesh.x-x).^2+(mesh.y-y).^2 )./obj.r0;
        end
        
        %> the exact function
        function f_ext = getExtFunc(obj, mesh, time)
            f_ext = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            
            temp = zeros( mesh.cell.Np, mesh.K );
            % The slotted cylinder
            theta = pi/2;
            [ x, y, r ]  = obj.getRelativeDistanceToExactPosition(mesh, theta, time);
            ind = ( r<=1.0) & ( (abs(mesh.x - x)>= 0.025) | (mesh.y-y>=0.1) );
            temp(ind) = 1;
            % The cone
            theta = -3*pi/4;
            [ ~, ~, r ]  = obj.getRelativeDistanceToExactPosition(mesh, theta, time);
            ind = ( r<=1.0);
            temp(ind) = 1 - r(ind);
            % The hump
            theta = -pi/4;
            [ ~, ~, r ]  = obj.getRelativeDistanceToExactPosition(mesh, theta, time);
            ind = ( r<=1.0);
            temp(ind) = (1 + cos( r(ind)*pi ))*0.25;
            
            f_ext(:,:,1) = temp;
            f_ext(:,:,2) = obj.w.* (- mesh.y);
            f_ext(:,:,3) = obj.w.*( mesh.x );
        end% func
    end
    
end

