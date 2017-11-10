function [ NCV, LAV ] = assembleLocalControlVolume( obj )

NCV = obj.Np;
LAV = zeros(NCV, 1);
r = obj.r;
s = obj.s;

for k = 1:obj.NFV
    vert = obj.EToV(:, k); % get the node index in each FV
    area = TriArea( r(vert), s(vert) );
    LAV(vert) = LAV(vert) + area/3;
end

end

