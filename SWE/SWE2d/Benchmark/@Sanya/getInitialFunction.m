function [ fext ] = getInitialFunction( obj )
fext = cell( obj.Nmesh, 1 );

for m = 1:obj.Nmesh%网格循环
    fext{m} = zeros(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield);%大小Np*K*Nfield
    topography_file =[pwd, '/SWE2d/@Sanya/mesh/bathymetry0111.txt'];
    fp = fopen(topography_file);
    fgets(fp);
    data = fscanf(fp, '%e %e %e', [3, inf]);
    fclose(fp);
    interp = scatteredInterpolant(data(1,:)',data(2,:)',...
        data(3,:)','linear');
    bot = - interp(obj.meshUnion(m).x, obj.meshUnion(m).y);
    h = max(0-bot,0);
    ind = ( any(h > 0) & any( bot > 0 ) );
    bot(:, ind) = 0 - h(:, ind);
    fext{m}(:,:,1) = h;
    fext{m}(:,:,4) = bot;
end

end%func
 


