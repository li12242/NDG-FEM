        function [x, rx, J] =GetEleGeometric(obj, vx)
            % get element Geometric Factor
            % Input:    vx  - Vertic Coordinate, size [2(nVertice) x nElement]
            % Output:   x   - node coordinate
            %           rx  - dr/dx at nodes
            %           J   - jacobi factor
            assert(size(vx,1)==2, 'transferToPhysic: input vx faults')
            x = 0.5*((1-obj.r)*vx(1,:) + (obj.r+1)*vx(2,:));
            xr  = obj.Dr*x; J = xr; rx = 1./J;
        end% function