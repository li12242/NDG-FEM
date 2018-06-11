%> \brief 2-dimensional non-linear shallow water equations
%> \details
%> This class descripe the conservation equations of the mass and
%> monuments of shallow water equations, written as
%> \f$ \frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot
%> \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U}), \f$
%> where \f$ \mathbf{U} = (h, hu, hv) \f$ are the conservative variables,
%> \f$ \mathbf{F}(\mathbf{U}) \f$ and \f$ \mathbf{S}(\mathbf{U}) \f$
%> are the flux terms and source term, respectively.
%> For the SWE model, the wet/dry (WD) probelm is addressed with the
%> methods fron Li (2018), which requires to determine the WD states of
%> each elements. The numerical flux
classdef SWEAbstract2d < NdgPhysMat
    
    properties (Abstract, Constant)
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
    end
    
    properties ( Constant )
        %> number of physical field
        Nfield = 5
        %> number of variable field
        Nvar = 3
        %> index of variable in physical field
        varFieldIndex = [ 1, 2, 3 ]
    end
    
    properties ( SetAccess = protected )
        %> gradient of bottom elevation
        zGrad
    end
    
    properties ( SetAccess = private )
        %> solver for coriolis source term
        coriolisSolver
        %> solver for friction source term
        frictionSolver
        %> solver for wind source term
        windSolver
        %> solver for unmerical flux
        numfluxSolver
        %> limiter type
        limiterSolver
    end
    
    % ======================================================================
    methods ( Hidden, Abstract ) % Abstract function, hidden
        %> abstract function to evaluate volume flux term
        [ E, G ] = matEvaluateFlux( obj, mesh, fphys );        
    end
    % ======================================================================

    
    % ======================================================================
    methods ( Abstract, Access = protected )
        %> determine wetting and drying status
        matUpdateWetDryState(obj, fphys)
        
        %> evaluate topography source term
        matEvaluateTopographySourceTerm( obj, fphys )
        
        %> evaluate post function
        [ fphys ] = matEvaluatePostFunc(obj, fphys)
    end
    % ======================================================================

    
    methods ( Hidden, Sealed, Access = public ) % public function, not allow to inherit
        
        %> impose boundary condition and evaluate cell boundary values
        [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext );
        [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext );
        
        %> evaluate local boundary flux
        function [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm )
            [ fluxM ] = mxEvaluateSurfFlux( obj.hmin, obj.gra, nx, ny, fm);
        end% func
        
        %> evaluate boundary numerical flux
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp )
            [ fluxS ] = obj.numfluxSolver.evaluate( obj.hmin, obj.gra, nx, ny, fm, fp );
        end% func
    end
        
    methods ( Sealed, Access = protected )
        [ fphys ] = matEvaluateLimiter( obj, fphys )
        
        %> determine time interval
        [ dt ] = matUpdateTimeInterval( obj, fphys )
        
        %> evaluate source term
        matEvaluateSourceTerm( obj, fphys )
    end
    
end
