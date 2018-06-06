function Fmask = assembleFacialNodeIndex(obj)
maxnfp = max(obj.Nfp);
Fmask = zeros(maxnfp, obj.Nface);
for f = 1:obj.Nface
    nfv = obj.Nfv(f);
    % get vertex index on face f
    rv = obj.vr( obj.FToV(1:nfv, f) );
    sv = obj.vs( obj.FToV(1:nfv, f) );
    tv = obj.vt( obj.FToV(1:nfv, f) );
    if(isrow(rv)) rv = rv'; end
    if(isrow(sv)) sv = sv'; end
    if(isrow(tv)) tv = tv'; end
    % get the nodes on face f
    cell = getStdCell(obj.N, obj.faceType(f));
    fr = cell.project_vert2node(rv);
    fs = cell.project_vert2node(sv);
    ft = cell.project_vert2node(tv);
    % get the nodes index
    for n = 1:obj.Nfp(f)
        dis = (fr(n) - obj.r).^2 + (fs(n) - obj.s).^2 + (ft(n) - obj.t).^2;
        ind = find(dis < 1e-10);
        Fmask(n, f) = ind;
    end
end
end% func