function [ stdCell ] = getStdCell( N, type )
%> @brief Factory of the StdCell subclasses.
%> @param[in] N    maximum order of the basis function
%> @param[in] type standard cell type, an enumStdCell member
%> @return stdCell the standard cell object
%> @note PrismTri is constructed with equal horizontal/vertical orders
%> (Nh = Nz = N); call StdPrismTri(Nh, Nz) directly for mixed orders.
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================
switch type
    case enumStdCell.Point
        stdCell = StdPoint(N);
    case enumStdCell.Line
        stdCell = StdLine(N);
    case enumStdCell.Tri
        stdCell = StdTri(N);
    case enumStdCell.Quad
        stdCell = StdQuad(N);
    case enumStdCell.PrismTri
        stdCell = StdPrismTri(N, N);
    otherwise
        error('getStdCell:unknownType', ...
            'Unsupported standard cell type: %s', char(type));
end
end% func
