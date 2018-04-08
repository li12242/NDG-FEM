classdef DamBreakSepUniformMesh2d < SWEPreBlanaced2d
    %DAMBREAKSEPUNIFORMMESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 1e-6
        %> gravity acceleration
        gra = 9.8
        %> Dam position
        damPosition = 0
        %> initial water depth
        hl = 5
        hr = 10
        ul = 0
        ur = 40
    end
    
    methods
        function obj = DamBreakSepUniformMesh2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
        end
        
        function verifySection( obj )
            Ng = 100;
            xg = linspace(-200, 400, Ng)'; yg = zeros(Ng, 1);
            pos = makeNdgPostProcessFromNdgPhys( obj );
            fphy = obj.fphys;
            fphyInterp = pos.interpolatePhysFieldToGaugePoint( fphy, xg, yg, yg );
            fext = obj.getExactFunction( obj.getOption('finalTime') );
            fextInterp = pos.interpolatePhysFieldToGaugePoint( fext, xg, yg, yg );
            figure('Color', 'w');
            subplot( 2, 2, [1,3] ); hold on; grid on; box on;
            plot( xg, fphyInterp(:,1), 'b.-' );
            plot( xg, fextInterp(:,1), 'r.-' );
            subplot( 2, 2, 2 ); hold on; grid on; box on;
            plot( xg, fphyInterp(:,2), 'b.-' );
            plot( xg, fextInterp(:,2), 'r.-' );
            subplot( 2, 2, 4 ); hold on; grid on; box on;
            uphyInterp = fphyInterp(:,2)./fphyInterp(:,1);
            uextInterp = fextInterp(:,2)./fextInterp(:,1);
            plot( xg, uphyInterp, 'b.-' );
            plot( xg, uextInterp, 'r.-' );
        end
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 6;
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
        
        function fphys = getExactFunction( obj, time )
            al = sqrt( obj.gra * obj.hl );
            cl = 2*al;
            Sl = obj.ul + cl;
            ar = sqrt( obj.gra * obj.hr );
            cr = 2*ar;
            Sr = obj.ur - cr;
            
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                atemp = zeros( mesh.cell.Np, mesh.K);
                utemp = zeros( mesh.cell.Np, mesh.K);
                mesh = obj.meshUnion(m);
                spe = (mesh.x - obj.damPosition)/time;
                
                ind = (spe <= obj.ul - al);
                atemp( ind ) = sqrt( obj.hl * obj.gra );
                utemp( ind ) = obj.ul;
                
                ind = ( (obj.ul - al) < spe ) & ( spe <= Sl );
                atemp( ind ) = ( obj.ul + cl - spe(ind) )/3;
                utemp( ind ) = ( obj.ul + cl + 2* spe(ind) )/3;
                
                ind = ( Sl <= spe ) & ( spe < Sr );
                atemp( ind ) = 0;
                utemp( ind ) = 0;
                
                ind = ( Sr <= spe ) & ( spe < (obj.ur + ar) );
                atemp( ind ) = ( spe(ind) - obj.ur + cr )/3;
                utemp( ind ) = ( 2*spe(ind) + obj.ur - cr )/3;
                
                ind = ( (obj.ur + ar) <= spe );
                atemp( ind ) = ar;
                utemp( ind ) = obj.ur;
                
                fphys{m}(:, :, 1) = atemp.^2 / obj.gra;
                fphys{m}(:, :, 2) = fphys{m}(:, :, 1) .* utemp;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.SlipWall, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-200, 400], [-10, 10], M, ceil(M/30), bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-200, 400], [-10, 10], M, ceil(M/30), bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

