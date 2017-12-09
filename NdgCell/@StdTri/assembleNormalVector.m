%> @brief Return the outword vector and the Jacobian factor of element boundary nodes.
%> @details
%> The outword normal vector transfer matrix is different with the
%> Jacobian transfer matrix. According to the differental chain
%> rule, the vector of \f$ (dx, dy) \f$ in physical domain can be
%> obtained from the equation
%> \f$ \begin{bmatrix} dx \cr dy \end{bmatrix} =
%> \mathbf{J} \cdot \begin{bmatrix} dr \cr ds \end{bmatrix} =
%> \begin{bmatrix} x_r & x_s \cr y_r & y_s \end{bmatrix} \cdot
%> \begin{bmatrix} dr \cr ds \end{bmatrix} \f$, where
%> \f$ (nr, ns) \f$ is the displacement in the reference element.
%>
%> However, as the vector \f$ (nx, ny) \f$ is normal to the
%> displacement vector \f$ \begin{bmatrix} nx \cr ny \end{bmatrix}^T
%> \cdot \begin{bmatrix} dx \cr dy \end{bmatrix} \f$,
%> the transfer matrix of the normal vector as
%> \f$ A = \left[ \mathbf{ J } ^{-1} \right]^T = \begin{bmatrix} y_s & -x_s \cr
%> -y_r & x_r \end{bmatrix}^T \Big/ \left| \mathbf{ J } \right|
%> = \begin{bmatrix} rx & ry \cr sx & sy \end{bmatrix} \f$.
%>
%> Therefore, the ourword vectors of the physical element are
%> derivaed from transforming the ourword vectors of the reference
%> element, while the Jacobian \f$ \mathbf{J}_s \f$ is obtained from the transforming
%> of edges vectors in the reference element.
function [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z )

nx = zeros( obj.TNfp, size(x, 2) );
ny = zeros( obj.TNfp, size(x, 2) );
nz = zeros( obj.TNfp, size(x, 2) );

vx = x( obj.Fmask(1,:), : );
vy = y( obj.Fmask(1,:), : );
% vz = z( obj.Fmask(1,:), : );

faceIndexStart = ones(obj.Nface, 1); % start index of each face node
for f = 2:obj.Nface
    faceIndexStart(f) = faceIndexStart(f-1) + obj.Nfp(f-1);
end
for f = 1:obj.Nface
    Nfp = obj.Nfp(f);
    face_x1 = vx( obj.FToV(1,f), : );
    face_x2 = vx( obj.FToV(2,f), : );
    face_y1 = vy( obj.FToV(1,f), : );
    face_y2 = vy( obj.FToV(2,f), : );
    
    ind = faceIndexStart(f):(faceIndexStart(f)+Nfp-1);
    nx(ind, :) = repmat( (face_y2 - face_y1), Nfp, 1 );
    ny(ind, :) = repmat(-(face_x2 - face_x1), Nfp, 1 );
end
Js = sqrt(nx.*nx+ny.*ny);
% normalise
nx = nx./Js;
ny = ny./Js;
Js = Js.*0.5;

% Nfp1 = 1:obj.Nfp(1);
% Nfp2 = (Nfp1(end)+1):sum( obj.Nfp(1:2) );
% Nfp3 = (Nfp2(end)+1):sum( obj.Nfp(1:3) );
% nr = zeros( obj.TNfp, 1 );
% ns = zeros( obj.TNfp, 1 );
% nr( Nfp1 ) = 0; ns( Nfp1 ) = -1;
% nr( Nfp2 ) = 1; ns( Nfp2 ) = 1;
% nr( Nfp3 ) = -1; ns( Nfp3 ) = 0;
%
% [ rx, ry, sx, sy, ~ ] = obj.assembleJacobianMatrix( x, y );
% nx = rx( obj.Fmask(:), : ) .* nr + sx( obj.Fmask(:), : ) .* ns;
% ny = ry( obj.Fmask(:), : ) .* nr + sy( obj.Fmask(:), : ) .* ns;
% len = sqrt( nx.^2 + ny.^2 );
% nx = nx ./ len;
% ny = ny ./ len;
% nz = zeros( size(nx) );

% xr = obj.Dr * x; xs = obj.Ds * x;
% yr = obj.Dr * y; ys = obj.Ds * y;
%
% Js = zeros( obj.TNfp, size( x, 2 ) );
% dr = 2; ds = 0;
% Fmask = obj.Fmask(:, 1);
% Js( Nfp1, : ) = sqrt( (dr*xr(Fmask, :) + ds*xs(Fmask, :) ).^2 + ( dr.*yr(Fmask, :) + ds*ys(Fmask, :) ).^2 );
% dr = -2; ds = 2;
% Fmask = obj.Fmask(:, 2);
% Js( Nfp2, : ) = sqrt( (dr*xr(Fmask, :) + ds*xs(Fmask, :) ).^2 + ( dr.*yr(Fmask, :) + ds*ys(Fmask, :) ).^2 );
% dr = 0; ds = -2;
% Fmask = obj.Fmask(:, 3);
% Js( Nfp3, : ) = sqrt( (dr*xr(Fmask, :) + ds*xs(Fmask, :) ).^2 + ( dr.*yr(Fmask, :) + ds*ys(Fmask, :) ).^2 );
% Js = Js ./ 2; % divide the length of standard line (=2)
end