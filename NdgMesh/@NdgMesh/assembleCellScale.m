function [LAV, charLength] = assembleCellScale( obj, J )

% Jacobian determination on each quadrature points
one = ones( obj.cell.Np, obj.K );
LAV = mxGetMeshIntegralValue( one, obj.cell.wq, J, obj.cell.Vq );
%LAV = mxGetIntegralValue(, obj.cell.wq, obj.Jq );
switch obj.cell.type
    case enumStdCell.Line
        charLength = LAV;
    case enumStdCell.Tri
        charLength = sqrt( 2*LAV );
    case enumStdCell.Quad
        charLength = sqrt( LAV );
end
end% func