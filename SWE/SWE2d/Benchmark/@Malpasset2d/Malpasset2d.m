classdef Malpasset2d < SWEPreBlanaced2d
    
    properties( Constant )
        hmin = 1e-1
        n = 0.029.^2
        gra = 9.81
        gmshfile = 'SWE/SWE2d/Benchmark/@Malpasset2d/mesh/malpasset.msh';
        trianglefile = 'mesh/triMeshLai/mesh';
    end
    
    properties
        %
        Ng
        % gauge points location
        xg, yg
        % 
        cellId
        %
        Vg
        %
        gaugeMaxDepthResult
    end
    
    methods(Access=protected, Hidden)
        function fphys = matInterpolateTopography( obj, fphys ) 
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                tpfile = ...
                     'SWE/SWE2d/Benchmark/@Malpasset2d/mesh/bathmetry1.txt';
                fp = fopen(tpfile);
                fgets(fp);
                data = fscanf(fp, '%e %e %e', [3, inf]);
                fclose(fp);
                interp = scatteredInterpolant( ...
                    data(1,:)', data(2,:)', data(3,:)', 'linear');

                fphys{m}(:,:,4) = interp(mesh.x, mesh.y);
            end
        end% func
    end
    
    methods
        function obj = Malpasset2d( N )
            obj = obj@SWEPreBlanaced2d();
            mesh = readTriMeshFile( N, obj.trianglefile );
            %mesh = makeGmshFileUMeshUnion2d( N, obj.gmshfile );
            obj.initPhysFromOptions( mesh );
            
            [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
            [ obj.cellId, obj.Vg ] = accessGaugePointMatrix( obj );
            obj.gaugeMaxDepthResult = zeros( obj.Ng, 1 );
        end
    end
    
    methods(Access=protected)
        function [ cellId, Vg ] = accessGaugePointMatrix( obj )
            cellId = cell( obj.Nmesh, 1 );
            Vg = cell( obj.Nmesh, 1 );
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                [ cellId{m}, Vg{m} ] = mesh.accessGaugePointLocation( ...
                    obj.xg, ...
                    obj.yg, ...
                    obj.xg );
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2000;
            
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('FrictionType') = FrictionType.Manning;
            option('FrictionCoefficient_n') = obj.n;
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            [ fphys ] = matEvaluatePostFunc@SWEAbstract2d( obj, fphys );
            for m = 1:obj.Nmesh
                for i = 1:obj.Ng
                    depth = obj.Vg{m}(i, :) * fphys{m}(:, obj.cellId{m}(i), 1);
                    obj.gaugeMaxDepthResult(i) = max( ...
                        obj.gaugeMaxDepthResult(i), ...
                        depth );
                end
            end
        end% func
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            end
                
            fphys = obj.matInterpolateTopography( fphys );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = fphys{m}(:, :, 4);
                
%                 h = 100 - fphys{m}(:, :, 4);
%                 ind = (mesh.EToR == int8(NdgRegionType.Dry) );
%                 h(:, ind) = 0;
                
                h = zeros(mesh.cell.Np, mesh.K);
                dam_x = [4701.18,4656.5]; 
                dam_y = [4143.41,4392.1];
                k = (dam_y(2)-dam_y(1))/(dam_x(2)-dam_x(1));
                flag = mesh.y - dam_y(1) - k*( mesh.x - dam_x(1) );
                I = find( flag <= 0 );
                h(I) = 100 - bot(I);
                
                I = ( abs( mesh.y - 5250 ) < 150 ) & ...
                    ( abs( mesh.x - 4500 ) < 200 );
                h(I) = 0;
%                 h(:, ind) = 0 - fphys{m}(:, ind, 4);
%                 h( h < 0 ) = 0;
                fphys{m}(:, :, 1) = h;
            end
        end% func
    end
    
end

function mesh = readTriMeshFile( N, casename )
[path, ~, ~] = fileparts( mfilename('fullpath') );
% vertex node file
nodefile = [path, '/', casename, '.node'];
fid1 = fopen(nodefile, 'r');

temp = fscanf(fid1,'%d %d %d %d',4);
Nv = temp(1);
data = fscanf(fid1,'%d %f %f\n',[3,Nv]);
vx = data(2, :)';
vy = data(3, :)';
fclose( fid1 );

% element file
elefile = [path, '/', casename, '.ele'];
fid1 = fopen(elefile, 'r');
temp = fscanf(fid1,'%d %d %d\n',3);
Ne = temp(1);
data = fscanf(fid1,'%d %d %d %d\n',[4,Ne]);
EToV = data(2:4,:);
fclose( fid1 );

EToR = ones( Ne, 1 );
tri = StdTri(N);
mesh = NdgMesh2d(tri, Nv, vx, vy, Ne, EToV, EToR, []);

temp = find(mesh.EToE == repmat( 1:mesh.K, mesh.cell.Nface, 1 ) );
[f, k] = ind2sub(...
    [mesh.cell.Nface, mesh.K], ...
    temp );
mesh.EToB( temp ) = NdgEdgeType.SlipWall;
for i = 1:numel(temp)
    finalNfp = sum( mesh.cell.Nfp(1:f(i)) );
    startNfp = finalNfp - mesh.cell.Nfp(f(i)) + 1;
    mesh.eidtype( startNfp:finalNfp, k(i) ) = int8( NdgEdgeType.SlipWall );
end

end% func

function [Ng, xg, yg] = setGaugePoint( )
xg = [5550, 11900, 13000, 4947.46, 5717.30, 6775.14,...
    7128.20, 8585.3, 9674.97, 10939.15, 11724.37, 12723.70...
    4913.11 5159.75 5790.63 5886.54 6763.05 6929.97 7326.02 7441.01...
    8735.94 8628.6 9761.13 9800 10957 11156.99 11689.05 11626.05 12333.72];

yg = [4400, 3250, 2700, 4289.71, 4407.61, 3869.23,...
    3162.00, 3443.08, 3085.89, 3044.78, 2810.41, 2485.08...
    4244.01 4369.62 4177.76 4503.97 3429.6 3591.87 2948.78 3232.12 3264.61...
    3604.63 3480.36 2414.79 2651.94 3800.72 2592.36 3406.8 2269.74];

Ng = numel(xg);
end
