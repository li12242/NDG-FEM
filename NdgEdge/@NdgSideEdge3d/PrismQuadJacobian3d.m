function [ nx, ny, nz, Js ] = PrismQuadJacobian3d( obj, mesh, f1, e1, fid )
    rx = mesh.rx( fid, e1 ); ry = mesh.ry( fid, e1 ); rz = mesh.rz( fid, e1 );
    sx = mesh.sx( fid, e1 ); sy = mesh.sy( fid, e1 ); sz = mesh.sz( fid, e1 );
    
    if f1 == 1
        nx = - sx;
        ny = - sy;
        nz = - sz;
    elseif f1 == 2
        nx = rx;
        ny = ry;
        nz = rz;
    elseif f1 == 3
        nx = sx;
        ny = sy;
        nz = sz;
    elseif f1 == 4
        nx = - rx;
        ny = - ry;
        nz = - rz;
    end
    
    Js = sqrt( nx .* nx + ny .* ny );
    nx = nx ./ Js; ny = ny ./ Js;
    Js = Js .* mesh.J( fid, e1 );
end