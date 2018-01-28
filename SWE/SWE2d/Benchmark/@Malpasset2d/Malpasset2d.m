classdef Malpasset2d < SWEPreBlanaced2d
    
    properties( Constant )
        hmin = 1e-1
        n = 0.029.^2
        gra = 9.81
        gmshfile = 'SWE/SWE2d/Benchmark/@Malpasset2d/mesh/malpasset.msh';
        trianglefile = 'mesh/triMeshLai/mesh';
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
        end
        
        function assessMexGaugeDepth( obj )
        end
    end
    
    methods(Access=protected)
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
        
%         function [ fphys ] = matEvaluatePostFunc(obj, fphys)
%             [ fphys ] = matEvaluatePostFunc@SWEAbstract2d( obj, fphys );
%             obj.matUpdateWetDryState( fphys );
%         end% func
%         
%         function fphys = matEvaluateLimiter( obj, fphys )            
%             obj.matUpdateWetDryState( fphys )
%             
%             fphys = obj.limiter.matLimit( fphys, 2 );
%             fphys = obj.limiter.matLimit( fphys, 3 );
%             for m = 1:obj.Nmesh % update new elevation
%                 fphys{m}(:,:,5) = fphys{m}(:,:,1) + fphys{m}(:,:,4);
%             end
%             fphys = obj.limiter.matLimit( fphys, 5 ); % enforce the elevation
%             for m = 1:obj.Nmesh % update new elevation
%                 fphys{m}(:,:,1) = fphys{m}(:,:,5) - fphys{m}(:,:,4);
%             end
%         end% func
        
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

temp = (mesh.EToE == repmat( 1:mesh.K, mesh.cell.Nface, 1 ) );
mesh.EToB( temp ) = int8(NdgEdgeType.SlipWall);
end
