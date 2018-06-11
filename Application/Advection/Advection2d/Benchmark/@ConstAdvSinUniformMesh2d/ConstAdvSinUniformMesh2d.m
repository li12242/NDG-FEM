classdef ConstAdvSinUniformMesh2d < AdvAbstractOBFlow2d
    %CONSTADVSINUNIFORMMESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> Number of basis function
        N
        %> Number of elements on each axis
        M
    end
    
    properties ( Constant )
        u0 = 0.5;
        v0 = 0;
        L = 1;
    end
    
    methods (Access = public)
        function obj = ConstAdvSinUniformMesh2d( N, M, cellType )
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@AdvAbstractOBFlow2d( );
            obj.N = N;
            obj.M = M;
            obj.initPhysFromOptions( mesh );
        end
    end

    methods( Access = protected )
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = getExtFunc(obj, mesh.x, mesh.y, 0);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 50;
            finalTime = 2.0;
            option('startTime') = 0.0;
            option('finalTime') = finalTime;
            option('timeInterval') ...
                = 2/obj.M/sqrt(obj.u0 ^ 2 + obj.v0 ^2)/(2*obj.N + 1);
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
        end
        
        function [ fext ] = getExtFunc( obj, x, y, time )
            fext = sin( 2 * pi * (x - obj.u0 * time) / obj.L );
        end
    end

    
end


function mesh = makeUniformMesh(N, M, type)
bctype = [...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], M, M, bctype);
else
    msgID = [mfilename, ':inputCellTypeError'];
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
