function [ visual ] = makeVisualizationFromNdgPhys( physMat )
    if ( physMat.meshUnion(1).type == enumMeshDim.One )
        [ visual ] = makeVisualizationFromNdgPhys1d( physMat );
    elseif ( physMat.meshUnion(1).type == enumMeshDim.Two )
        [ visual ] = makeVisualizationFromNdgPhys2d( physMat );
    end
end

function [ visual ] = makeVisualizationFromNdgPhys1d( physMat )
    visual = [];
    for m = 1 : physMat.Nmesh
        visual = [ visual, Visual2d( physMat.meshUnion(m) ) ];
    end
end

function [ visual ] = makeVisualizationFromNdgPhys2d( physMat )
    visual = [];
    for m = 1 : physMat.Nmesh
        visual = [ visual, Visual2d( physMat.meshUnion(m) ) ];
    end
end