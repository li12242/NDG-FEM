%> @brief Abstract class of the solver.
%
%> NdgPhys is the abstract superclass of all the solvers. The class stores
%> the options for solving problem with specific spacial scheme, such as 
%> using the quadrature-free scheme for strong form discrete equations, 
%> and the temporal discrete methods used.
%>
%> The public functions includes: 
%> @code
%>  [ ] = initPhysFromOptions( obj, mesh ); // Public function to call initial function;
%>  [ ] = listOption( obj ); // List all the options in the solver;
%>  [ field ] = getOption( obj, fieldName ); // get the option values;
%>  [ ] = addOption( obj, optionName, optionValue ); // add new options into the solver; 
%>  [ ] = draw(obj, fieldId); // draw the field;
%> @endcode
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
        %> CFL number
        cfl
    end
    
    properties( SetAccess = protected )
        %> setting options
        option = containers.Map();
    end
    
    methods( Access = public )
        
        %> @brief Initial function
        %> This function will set the solver's mesh objects and initialize
        %> the physical field
        initPhysFromOptions( obj, mesh )

        %> @brief Add new setting option into the solver
        function addOption( obj, optionName, optionValue )
            obj.option( optionName ) = optionValue;
        end
        
        %> @brief List all the options in the solver
        function listOption( obj )
            keys = obj.option.keys;
            values = obj.option.values;
            fprintf('Options = {\n');
            for n = 1:obj.option.Count
                if isnumeric( values{n} )
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
                msgID = [ mfilename, ':UnknownOptionField'];
                msgtext = ['The option field ', fieldName ,' does not exist.'];
                throw( MException(msgID, msgtext) );
            end
        end% func
        
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
    
end
