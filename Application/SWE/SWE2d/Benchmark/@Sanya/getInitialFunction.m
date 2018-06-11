function [ fext ] = getInitialFunction( obj )

fext = cell( obj.Nmesh, 1 );

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    fext{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
    topography_file = ...
        [ fileparts( mfilename('fullpath') ) '/mesh/bathymetry0111.txt' ];
    data = load( topography_file );

    interp = scatteredInterpolant(data(:, 1),data(:, 2), data(:, 3),'linear');
    bot = - interp( mesh.x, mesh.y);
    h = max( 0 - bot, 0 );
%     ind = ( any(h > 0) & any( bot > 0 ) );
%     bot(:, ind) = 0 - h(:, ind);
    fext{m}(:,:,1) = h;
    fext{m}(:,:,4) = bot;
end

end%func
 


