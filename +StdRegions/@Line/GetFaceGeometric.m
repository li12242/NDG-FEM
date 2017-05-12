  function [nx, sJ] = GetFaceGeometric(~, x)
            % get Face Normal vector & surface jacobi factor
            % Input:    x  - node coordinate, size [nNode, nElement]
            % Output:   nx - outward normal vector
            %           sJ - surface edge jacobi coefficient
            nx = zeros(2, size(x,2));
            % Define outward normals
            nx(1, :) = -1.0; nx(2, :) = 1.0; sJ = ones(size(nx));
        end% function