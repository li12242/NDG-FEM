function [ mesh ] = makeGmshFileUMeshUnion2d( N, filename )

fid1 = fopen(filename, 'r');
if( fid1 < 0 )
    msgID = [mfilename, ':inputFileNameError'];
    msgtext = ['The input file name: ', filename...
        ' is incorrect'];
    ME = MException(msgID, msgtext);
    throw(ME);
end

for i=1:4
   fgetl(fid1);
end
Nv = fscanf(fid1,'%d',1);
data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
vx = data(2, :)';
vy = data(3, :)';
for i=1:2
   fgetl(fid1);
end
Ne = fscanf(fid1,'%d\n',1);
Nedge = 0; 
Ntri = 0; 
Nquad = 0;
for i = 1:Ne
    temp = fgetl(fid1);
    data = str2num(temp);
    switch data(2)
        case 1
            Nedge = Nedge + 1;
        case 2
            Ntri = Ntri + 1;
        case 3
            Nquad = Nquad + 1;
    end
end

fseek(fid1,0,-1);
for i=1:(Nv+8)
    fgetl(fid1);
end

data = fscanf(fid1,'%d %d %d %d %d %d %d\n',[7,Nedge]);
BCToV = [data(6, :); data(7, :); data(4, :)];

if Ntri > 0
    data = fscanf(fid1,'%d %d %d %d %d %d %d %d \n',[8, Ntri]);
    EToVTri = [data(6, :); data(7, :); data(8, :)];
    EToRTri = data(4, :);
    stdTri = StdTri(N);
    triMesh = NdgMesh2d(stdTri, Nv, vx, vy, Ntri, EToVTri, EToRTri, BCToV);
end

if Nquad > 0
    data = fscanf(fid1,'%d %d %d %d %d %d %d %d %d\n',[9, Nquad]);
    EToVQuad = [data(6, :); data(7, :); data(8, :); data(9, :)];
    EToRQuad = data(4, :);
    stdQuad = StdQuad(N);
    quadMesh = NdgMesh2d(stdQuad, Nv, vx, vy, Nquad, EToVQuad, EToRQuad, BCToV);
end

if (Ntri > 0) && (Nquad > 0)
    mesh = makeMeshUnion(2, triMesh, quadMesh);
elseif (Ntri > 0) && (Nquad == 0)
    mesh = triMesh;
elseif (Nquad > 0) && (Ntri == 0)
    mesh = quadMesh;
end

end

