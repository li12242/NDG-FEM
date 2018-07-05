classdef LSWEOutput3d < NcOutput
    %LSWEOUTPUT3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant )
        Nfield3d = 7;
        Nfield2d = 5;
    end
    
    properties ( SetAccess = protected )
        mesh2d
        mesh3d
        fieldId2
        fieldId3
    end
    
    methods
        function obj = LSWEOutput3d( casename, Nfield, dt )
            obj = obj@NcOutput( casename, Nfield, dt );
        end
        
        function initFromMesh( obj, mesh2d, mesh3d )
            % define dimension
            dimTime = NdgNcDim('Nt', 0);
            dimK2 = NdgNcDim('K2d', mesh2d.K);
            dimNp2 = NdgNcDim('Np2d', mesh2d.cell.Np);
            dimNfield2 = NdgNcDim('Nfield2d', obj.Nfield2d );
            
            dimK3 = NdgNcDim('K3d', mesh3d.K );
            dimNp3 = NdgNcDim('Np3d', mesh3d.cell.Np );
            dimNfield3 = NdgNcDim('Nfield3d', obj.Nfield3d );
            
            % define variable
            varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
            varField2 = NdgNcVar('fphys2d', ...
                [ dimNp2, dimK2, dimNfield2, dimTime], ...
                enumNcData.NC_DOUBLE);
            
            varField3 = NdgNcVar('fphys3d', ...
                [ dimNp3, dimK3, dimNfield3, dimTime], ...
                enumNcData.NC_DOUBLE);
            
            obj.filename = [ obj.casename, '/', obj.casename, '.nc' ];
            obj.ncfile = NdgNcFile( obj.filename, ...
                [dimTime, dimK2, dimNp2, dimNfield2, dimK3, dimNp3, dimNfield3], ...
                [varTime, varField2, varField3]);
            
            if ~isdir(obj.casename)
                mkdir(obj.casename);
            end
            % init file
            obj.ncfile.defineIntoNetcdfFile();

            % set properties
            obj.timeVarableId = varTime.id;
            obj.fieldId2 = varField2.id;
            obj.fieldId3 = varField3.id;
            obj.mesh2d = mesh2d;
            obj.mesh3d = mesh3d;
        end
        
        function outputResult( obj, time, field2d, field3d )
            if ( time - obj.timePrevious ) > obj.timeInterval
                % output time
                startInd = obj.outputStep;
                countInd = 1;
                netcdf.putVar(obj.ncfile.ncid, obj.timeVarableId, startInd, countInd, time);

                % output physical field
                startInd = [ 0, 0, 0, obj.outputStep ];
                countInd = [ obj.mesh2d.cell.Np, obj.mesh2d.K, obj.Nfield2d, 1 ];
                netcdf.putVar(obj.ncfile.ncid, obj.fieldId2, startInd, countInd, field2d);

                startInd = [ 0, 0, 0, obj.outputStep ];
                countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, obj.Nfield3d, 1 ];
                netcdf.putVar(obj.ncfile.ncid, obj.fieldId3, startInd, countInd, field3d);

                % increase output step num
                obj.outputStep = obj.outputStep + 1;
                obj.timePrevious = time;
            end
        end
        
        function [fphys2d, fphys3d] = readOutputResult( obj, timeStep )
            if timeStep <= obj.outputStep                
                startInd = [ 1, 1, 1, timeStep ];
                countInd = [ obj.mesh2d.cell.Np, obj.mesh2d.K, obj.Nfield2d, 1 ];
                fphys2d = ncread( obj.filename, 'fphys2d', startInd, countInd );
                
                startInd = [ 1, 1, 1, timeStep ];
                countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, obj.Nfield3d, 1 ];
                fphys3d = ncread(obj.filename, 'fphys3d', startInd, countInd );
            else
                error('Output time step is less than inptu step number!')
            end
        end
    end
    
end

