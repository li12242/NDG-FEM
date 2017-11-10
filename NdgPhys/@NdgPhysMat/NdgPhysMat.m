%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgPhysMat < NdgPhys
    
    properties( SetAccess = protected )
        %> cell array for RHS
        frhs
        %> cell array for external values
        fext
        %> output netcdf file objects
        outputNcFile
        %> limiter object
        limiter
    end
    
    properties( SetAccess = protected )
        ftime
        outputStep
        outputStepInterval
        outputTimeInterval
    end
    
    methods
        function obj = NdgPhysMat()
            obj = obj@NdgPhys();
        end
        
        %> solve with the Matlab function
        function matSolve( obj )
            obj.matEvaluateTemporalDiscrete();
        end
        
        %> @brief Public function to call initial function
        function setNdgPhys( obj, mesh )
            setNdgPhys@NdgPhys( obj, mesh );
            for m = 1:obj.Nmesh
                Np = obj.meshUnion(m).cell.Np;
                K = obj.meshUnion(m).K;
                obj.frhs{m} = zeros( Np, K, obj.Nvar );
            end
            [ obj.fext ] = obj.setInitialField( );
        end% func
    end
    methods( Abstract, Access = protected )
        %> @brief A warper functions for evaluating the RHS term
        %> @details
        %> For the 2d problem, the function should call function
        %> matEvaluateRHS2d; For the 3d problem, the function
        %> matEvaluateRHS3d should be called.
        matEvaluateRHS( obj, fphys )
        %> @brief function for calculating the flux term
        %> @details
        %> Function to calculate the flux term of the equation.
        %> @param[in] mesh The ith mesh object
        %> @param[in] fphys The physical field on the ith mesh grid
        %> @retval[out] the flux term of the variable on each axis
        matEvaluateFlux( obj, mesh, fphys )
        %> @brief function for calculating the flux division
        %> @details
        %> Function calculates the flux division of the normal flux
        %> and the numerical flux term on mesh, which is obtained from
        %> \f$ \mathbf{F} \cdot \mathbf{n} - F^*  \f$
        %> @param [in] mesh The ith mesh object
        %> @param [in] fphys The physical field on the mesh grid
        %>
        matEvaluateNumericalFlux( obj, mesh, fphys, fext )        
    end
    
    methods( Access = protected )
                
        %> @brief Temporal discrete function
        %> @details
        %> The temporal discrete function will solve the function with the
        %> specific temporal discrete methods, which is defined in option
        %> by the name 'temporalDiscreteType'.
        function matEvaluateTemporalDiscrete( obj )
            obj.matInitMatSolver();
            
            switch obj.getOption('temporalDiscreteType')
                case NdgTemporalDiscreteType.Euler
                    % call the Euler temporal discrete function
                    obj.matEvaluateEuler();
                case NdgTemporalDiscreteType.RK45
                    % call the SSP-RK45 temporal discrete function
                    obj.matEvaluateRK45();
                otherwise
                    msgID = 'NdgPhysMat:UnknownTemproalDicsreteType';
                    msgtext = ['The temporal discrete type ', ...
                        obj.getOption('temporalDiscreteType') ,' is invalid.'];
                    throw( MException(msgID, msgtext) );
            end
        end% func
        
        %> An instantiation of the temporal-discrete function with the SSP-RK45 method
        matEvaluateRK45( obj );
        %> An instantiation of the temporal-discrete with the Euler method
        matEvaluateEuler( obj );
        
        
        matUpdateExternalField( obj, time, fphys )
        dt = matUpdateTimeInterval( obj, fphys )
        matEvaluateRHS2d( obj, fphys )
        matEvaluateRHS3d( obj, fphys )
        matEvaluateSourceTerm( obj, fphys )
        fphys = matEvaluateLimiter( obj, fphys )
        fphys = matEvaluatePostFunc( obj, fphys )
        
        %> @brief
        matInitMatSolver( obj )
        matUpdateOutputResult( obj, time, step, fphys )
        matCloseOutputNetcdfFile( obj )
    end
    
end
