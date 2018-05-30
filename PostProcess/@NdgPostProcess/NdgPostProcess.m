classdef NdgPostProcess < handle
    
    properties
        %> number of mesh
        Nmesh
        %> mesh objects
        meshUnion
        %> output NetCDF file name
        outputFile
        %> number of variable field
        Nvar
        %> number of output step
        Nt
        %> output time
        time
    end
    
    methods
        function obj = NdgPostProcess( meshUnion, casename )
            [ obj.Nmesh ] = numel( meshUnion );
            [ obj.meshUnion ] = meshUnion;
            [ obj.outputFile ] = cell( obj.Nmesh, 1 );
            
            for n = 1:obj.Nmesh
                [ obj.outputFile{n} ] = [ casename, '.', num2str(n), '-', ...
                    num2str(obj.Nmesh), '.nc'];
            end
            
            [ obj.Nt ] = accessOutputStepNumber( obj );
            [ obj.Nvar ] = accessOutputVarNumber( obj );
            [ obj.time ] = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                obj.time{m} = obj.assessOutputVar( m, 'time' );
            end
        end
        
        %======================================================================
        %> \brief draw output physical field
        %>
        %> More detailed description.
        %>
        %> \param physField external field
        %======================================================================
        function drawResult( obj, delta, fieldId, physField )
            varId = fieldId;
            for t = 1:delta:obj.Nt
                field = obj.accessOutputResultAtStepNum( t );
                for m = 1:obj.Nmesh
                    obj.meshUnion(m).draw( field{m}(:,:,varId) + physField );
                end
                drawnow;
            end
        end
        
        [ mass ] = checkMassVolume( obj, varId )
        [ err ] = evaluateNormErr1( obj, fphys, fext );
        [ err ] = evaluateNormErr2( obj, fphys, fext );
        [ err ] = evaluateNormErrInf( obj, fphys, fext );
        
        [ fg ] = interpolateOutputStepResultToGaugePoint( obj, xg, yg, zg, outputStep );
        [ fg ] = interpolateOutputResultToGaugePoint( obj, xg, yg, zg );
        [ fg ] = interpolatePhysFieldToGaugePoint( obj, fphys, xg, yg, zg );
        
        [ Noutput ] = accessOutputStepNumber( obj )
        [ fphys ] = accessOutputResultAtStepNum( obj, stepId )
        [ Nvar ] = accessOutputVarNumber( obj )
        
        function [ var ] = assessOutputVar( obj, meshId, varName )
            var = ncread(obj.outputFile{meshId}, varName);
        end
    end
    
end

