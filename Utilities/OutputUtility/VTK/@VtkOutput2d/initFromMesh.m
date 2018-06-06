function initFromMesh( obj, mesh )
%initFromMesh - Description
%
% Syntax: output = initFromMesh(obj)
%
% Long description

    [ obj.Np, Ncell, Ncon, obj.SEToV, ctype ] = InitSubConnect( mesh );

    obj.CellVertList = zeros( obj.Np, Ncell * mesh.K );

    sk = 1;
    for k = 1 : mesh.K
        nodesk = mesh.cell.Np * ( k - 1 );
        for m = 1 : Ncell
            obj.CellVertList( :, sk ) = nodesk + obj.SEToV(:, m) - 1; % transform to 0-offset
            sk = sk + 1;
        end
    end

    obj.Npoint = mesh.K * mesh.cell.Np;
    obj.Points = [ mesh.x(:), mesh.y(:), zeros(obj.Npoint, 1) ]';
    obj.Ncell = mesh.K * Ncell;
    obj.Ncon = mesh.K * Ncon * Ncell;
    obj.ctype = repmat( int8(ctype), obj.Ncell, 1 );
    
    if ~isdir(obj.casename)
        mkdir(obj.casename);
    end
end

function [ Np, Ncell, Ncon, EToV, ctype ] = InitSubConnect( mesh )
    if ( mesh.cell.type == enumStdCell.Tri )
        [ Ncell, EToV ] = InitTriConnect2d( mesh.cell.N );
        Np = 3;
        Ncon = 3 + 1;
        ctype = enumVtkCell.VTK_TRIANGLE;
    elseif ( mesh.cell.type == enumStdCell.Quad )
        [ Ncell, EToV ] = InitQuadConnect2d( mesh.cell.N );
        Np = 4;
        Ncon = 4 + 1;
        ctype = enumVtkCell.VTK_QUAD;
    end
end