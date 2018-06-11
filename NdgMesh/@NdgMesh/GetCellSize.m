function obj = GetCellSize( obj )

% Jacobian determination on each quadrature points
one = ones( obj.cell.Np, obj.K );
LAV = mxGetMeshIntegralValue( one, obj.cell.wq, obj.J, obj.cell.Vq );

switch obj.cell.type
    case enumStdCell.Line
        obj.charLength = LAV;
    case enumStdCell.Tri
        obj.charLength = sqrt( 2*LAV );
    case enumStdCell.Quad
        obj.charLength = sqrt( LAV );
end

obj.LAV = LAV;

end% func