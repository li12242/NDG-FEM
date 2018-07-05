function obj = ExtendMesh3d( obj, cell, mesh2d, Nz )
    Kh = mesh2d.K;
    Kloc = obj.K;

    % extend element-vertex connection
    obj.EToV = [];
    for n = 1 : Nz
        obj.EToV = [ obj.EToV; 
            mesh2d.EToV + mesh2d.Nv * n;
            mesh2d.EToV + mesh2d.Nv * (n - 1) ];
    end
    obj.EToV = reshape( obj.EToV, cell.Nv, Kloc );

    % extend element connection
    obj.EToE = [];
    for n = 1 : Nz
        ind = ( (1 : Kh) - 1) * Nz + n;
        obj.EToE = [ obj.EToE; (mesh2d.EToE - 1) * Nz + n ];
        
        if n == Nz % bottom layer
            obj.EToE = [ obj.EToE; ind ];
        else
            obj.EToE = [ obj.EToE; ind + 1 ];
        end
        
        if n == 1 % top layer
            obj.EToE = [ obj.EToE; ind ];
        else
            obj.EToE = [ obj.EToE; ind - 1 ];
        end

    end% for
    obj.EToE = reshape( obj.EToE, cell.Nface, Kloc );

    obj.EToF = [];
    for n = 1 : Nz
        topId = ones( 1, mesh2d.K ) * 5;
        botId = ones( 1, mesh2d.K ) * 4;
        obj.EToF = [ obj.EToF; mesh2d.EToF ];
        
        if n == Nz % bottom layer
            obj.EToF = [ obj.EToF; botId ];
        else
            obj.EToF = [ obj.EToF; topId ];
        end
        
        if n == 1 % top layer
            obj.EToF = [ obj.EToF; topId ];
        else
            obj.EToF = [ obj.EToF; botId ];
        end

    end% for
    obj.EToF = reshape( obj.EToF, cell.Nface, Kloc );

    % extend mesh
    obj.EToM = [];
    ind = ones( 1, mesh2d.K ) * obj.ind;
    for n = 1 : Nz
        obj.EToM = [ obj.EToM; mesh2d.EToM; ind; ind ];
    end% for
    obj.EToM = reshape( obj.EToM, cell.Nface, Kloc );
    
    obj.EToL = zeros( Kh, Nz );
    for n = 1 : Nz
        obj.EToL(:, n) = ( (1 : Kh) - 1) * Nz + n;
    end
end