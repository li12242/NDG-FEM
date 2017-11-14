%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgPhys < handle
    
    properties( SetAccess = protected )
        %> cell array for physical field variable
        fphys
    end
    
    properties( Abstract, Constant )
        %> number of physical field
        Nfield
        %> number of variable field
        Nvar
        %> index of variable in physical field
        varFieldIndex
    end
    
    properties( SetAccess = protected )
        %> mesh objects
        meshUnion
        %> number of mesh
        Nmesh
    end
    
    properties( SetAccess = protected )
        %> setting options
        option
    end
    
    methods
        function obj = NdgPhys()
        end
        
        %> @brief Public function to call initial function
        function initPhysFromOptions( obj, mesh )
            [ obj.meshUnion ] = mesh;
            [ obj.Nmesh ] = numel(mesh);
            [ obj.option ] = obj.setOption( containers.Map() );
            [ obj.fphys ] = obj.setInitialField( );
        end% func
        
        %> @brief List all the options in the solver
        function listOption( obj )
            keys = obj.option.keys;
            values = obj.option.values;
            fprintf('option = {\n');
            for n = 1:obj.option.Count
                if isnumeric(values{n})
                    fprintf('\t %s: %d\n', keys{n}, values{n});
                else
                    fprintf('\t %s: %s\n', keys{n}, values{n});
                end
            end% for
            fprintf('}\n');
        end% func
        
        %> Get option from the option properties.
        function value = getOption( obj, fieldName )
            if obj.option.isKey( fieldName )
                value = obj.option( fieldName );
            else
                msgID = 'NdgPhys:UnknownOptionField';
                msgtext = ['The option field ', fieldName ,' is invalid.'];
                throw( MException(msgID, msgtext) );
            end
        end% func
        
%         [ err ] = evaluateNormErr2( obj );
%         [ err ] = evaluateNormErr1( obj );
%         [ err ] = evaluateNromErrInf( obj );
%         
%         [ fg ] = interpolateOutputStepResultToGaugePoint( obj, xg, yg, zg, outputStep );
%         [ fg ] = interpolateOutputResultToGaugePoint( obj, xg, yg, zg );
%         [ fg ] = interpolatePhysFieldToGaugePoint( obj, xg, yg, zg );
        
        %> @brief Draw the physical field on all meshes.
        function draw(obj, fieldId)
            for m = 1:obj.Nmesh
                obj.meshUnion(m).draw( obj.fphys{m}(:,:,fieldId) );
            end
        end% func
    end
    
    methods( Abstract, Access = protected )
        %> Set the solver options
        [ option ] = setOption( obj, option );
        %> Set the initial field
        [ fphys ] = setInitialField( obj );
    end
    
    methods( Access = protected )
%         [ Noutput ] = accessOutputResultStepNumber( obj )
%         [ fphys ] = accessOutputResultAtStepNum(obj, stepId)
    end
end
