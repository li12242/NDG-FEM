function [ adv, vis ] = initSolver( physMat )
%INITSOLVER Summary of this function goes here
%   Detailed explanation goes here

integralType = physMat.getOption('integralType');
equType = physMat.getOption('equationType');
dim = physMat.meshUnion(1).type;

if ( dim == enumMeshDim.One )
    adv = initAdvSolver1d( physMat, integralType, equType );
elseif ( dim == enumMeshDim.Two )
    adv = initAdvSolver2d( physMat, integralType, equType );
elseif ( dim == enumMeshDim.Three )
    adv = initAdvSolver3d( physMat, integralType, equType );
end

vis = NdgNonVisSolver( physMat );
end

function [ adv ] = initAdvSolver3d( physMat, integralType, equType )

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgQuadFreeStrongFormAdvSolver3d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgQuadFreeWeakFormAdvSolver3d( physMat );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgGaussQuadStrongFormAdvSolver3d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgGaussQuadWeakFormAdvSolver3d( physMat );
    end
end

end

function [ adv ] = initAdvSolver1d( physMat, integralType, equType )

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgQuadFreeStrongFormAdvSolver1d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgQuadFreeWeakFormAdvSolver1d( physMat );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgGaussQuadStrongFormAdvSolver1d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgGaussQuadWeakFormAdvSolver1d( physMat );
    end
end

end

function [ adv ] = initAdvSolver2d( physMat, integralType, equType )

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgQuadFreeStrongFormAdvSolver2d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgQuadFreeWeakFormAdvSolver2d( physMat );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgGaussQuadStrongFormAdvSolver2d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgGaussQuadWeakFormAdvSolver2d( physMat );
    end
end

end


