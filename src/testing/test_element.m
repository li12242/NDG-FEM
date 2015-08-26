function test_element

for degree = 1:5
    shape = femlib.element_type('dim', 1, 'vertice', 2, 'degree', degree, 'name', 'line');
    compare_line(shape)
end

% triangle test
% face = face_element_type('dim', 2, 'vertice', 3, 'degree', 3, 'name', 'triangle');
for degree = 1:5
    shape = femlib.element_type('dim', 2, 'vertice', 3, 'degree', degree, 'name', 'triangle');
    compare_tri(shape)
end

    function compare_line(shape)
        N = shape.nDegree;
        Np = N+1; Nfp = 1; Nfaces=2;
        r = JacobiGL(0,0,N);
        V  = Vandermonde1D(N, r); invV = inv(V);
        MassMatrix = invV'*invV;
        compare(shape.MassMatrix, MassMatrix, ['Line element MassMatrix test, Degree ',num2str(N)])
        % compare DifferenceMatrix
        Dr = Dmatrix1D(N, r, V);
        compare(shape.DifferenceMatrix, Dr, ['Line element DifferenceMatrix test, Degree ',num2str(N)])
        % compare EdgeMatrix
        Emat = zeros(Np,Nfaces*Nfp);
        Emat(1,1) = 1.0; Emat(Np,2) = 1.0;
        compare(shape.EdgeMassMatrix, Emat, ['Line element EdgeMassMatrix test, Degree ',num2str(N)])
    end %compare_line

    function compare_tri(shape)
        N = shape.nDegree;
        Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;
        [x,y] = Nodes2D(N); [r,s] = xytors(x,y);
        V = Vandermonde2D(N,r,s); invV = inv(V);
        MassMatrix = invV'*invV;
        compare(shape.MassMatrix, MassMatrix, ['Triangle element MassMatrix test, Degree ',num2str(N)])
        % compare DifferenceMatrix
        [Dr,Ds] = Dmatrices2D(N, r, s, V);
        DifferenceMatrix(:,:,1) = Dr; DifferenceMatrix(:,:,2) = Ds; 
        compare(shape.DifferenceMatrix, DifferenceMatrix, ['Triangle element DifferenceMatrix test, Degree ',num2str(N)])
        
        % compare EdgeMatrix
        Emat = zeros(Np, Nfaces*Nfp);
        fmask1   = find( abs(s+1) < NODETOL)'; 
        fmask2   = find( abs(r+s) < NODETOL)';
        fmask3   = find( abs(r+1) < NODETOL)';
        Fmask  = [fmask1;fmask2;fmask3]';
        % face 1
        faceR = r(Fmask(:,1));
        V1D = Vandermonde1D(N, faceR); 
        massEdge1 = inv(V1D*V1D');
        Emat(Fmask(:,1),1:Nfp) = massEdge1;
        % face 2
        faceR = r(Fmask(:,2));
        V1D = Vandermonde1D(N, faceR);
        massEdge2 = inv(V1D*V1D');
        Emat(Fmask(:,2),Nfp+1:2*Nfp) = massEdge2;
        % face 3
        faceS = s(Fmask(:,3));
        V1D = Vandermonde1D(N, faceS); 
        massEdge3 = inv(V1D*V1D');
        Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;
        compare(shape.EdgeMassMatrix, Emat, ['Triangle element EdgeMassMatrix test, Degree ',num2str(N)])
    end %function compare

    function compare(a,b, testName)
        TOL = 1e-12;    % global tol
        ndim = ndims(a);
        switch ndim
            case 1
                bool = all(abs(a - b) < TOL);
            case 2 
                bool = all(all(abs(a - b) < TOL));
            case 3
                bool = all(all(all(abs(a - b) < TOL)));
        end
        if bool
            fprintf(['pass: [test:',testName,']\n'])
        else
            error(['fail: [test:',testName,']\n'])
        end
    end %function compare
end