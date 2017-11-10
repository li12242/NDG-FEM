function [ EToE, EToF, EToM ] = assembleCellConnect(obj)
Nface = obj.cell.Nface;
Ke = obj.K;
faceId = obj.assembleGlobalFaceIndex();
EToE = ones(Nface, 1)*(1:Ke);
EToF = (1:Nface)'*ones(1,Ke);

for n = 1:( Nface*Ke )
    m = find( abs(faceId - faceId(n))<1e-10 );
    t = m( m ~= n );
    if( ~isempty(t) )
        [f, k] = ind2sub([ Nface, Ke ], t);
        EToE(n) = k;
        EToF(n) = f;
        %etoe(n) = fix( (t-1)./obj.cell.Nface )+1;
        %etof(n) = rem(t-1, obj.cell.Nface)+1;
    end
end

EToM = ones(Nface, Ke);
end