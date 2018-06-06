function obj = InitVisual2d(obj)

    % init SEToV
    if ( obj.mesh.cell.type == enumStdCell.Tri )
        [ obj.Ncell, obj.SEToV ] = InitTriVisual2d( obj.mesh.cell.N );
    elseif ( obj.mesh.cell.type == enumStdCell.Quad )
        [ obj.Ncell, obj.SEToV ] = InitQuadVisual2d( obj.mesh.cell.N );
    end

    obj.Ntri = obj.Ncell * obj.mesh.K;

    obj.tri = zeros( obj.Ntri, 3 );
    sk = 1;
    for i = 1 : obj.mesh.K
        vert = (i - 1) * obj.mesh.cell.Np;
        for j = 1 : obj.Ncell
            obj.tri(sk, :) = vert + obj.SEToV(:, j);
            sk = sk + 1;
        end
    end
end

function [ Ncell, EToV ] = InitTriVisual2d( N )
    Np = ( N + 1 ) * ( N + 2 ) / 2;
    Ncell = N^2;
    EToV = zeros(3, Ncell);

    sk = 1; 
    s2 = Np; 
    s1 = Np - 2;
    for row = 1 : N
        v1 = s2; v2 = s1;
        for kb = 1:row
            EToV(:, sk) = [v1, v2, v2+1]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end
        
        v1 = s2; v2 = s1;
        for kt = 1:row-1
            EToV(:, sk) = [v1, v2+1, v1+1]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end

        s1 = s1 - (row + 2);
        s2 = s2 - (row + 1);
    end

end

function [ NFV, EToV ] = InitQuadVisual2d( N )
    Np = ( N + 1 ) .^ 2;
    NFV = N ^ 2 * 2;
    EToV = zeros(3, NFV);

    sk = 1; 
    s2 = Np - N; 
    s1 = s2 - ( N + 1 );
    for row = 1 : N
        v1 = s1; v2 = s2;
        for kb = 1:N
            EToV(:, sk) = [ v1, v1 + 1, v2 ]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end
        
        v1 = s1; v2 = s2;
        for kt = 1 : N
            EToV(:, sk) = [ v1 + 1, v2 + 1, v2 ]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end

        s1 = s1 - (N + 1);
        s2 = s2 - (N + 1);
    end

end