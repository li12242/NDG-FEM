%> @brief Evaluating the RHS term for the 2d problem
%> @details The function calculate the RHS term on each mesh
function matEvaluateRHS2d( obj, fphys )

for m = 1:obj.Nmesh % calculate RHS term on each mesh
    mesh = obj.meshUnion(m);
    [ E, G ] = obj.matEvaluateFlux( mesh, fphys{m} );
    [ dFlux ] = obj.matEvaluateNumericalFlux( mesh, fphys{m}, obj.fext{m} );
    
    for i = 1:obj.Nvar
        %i = obj.varFieldIndex(i);
        [ obj.frhs{m}(:,:,i) ] = ...
            - mesh.rx.*( mesh.cell.Dr * E(:,:,i) ) ...
            - mesh.sx.*( mesh.cell.Ds * E(:,:,i) ) ...
            - mesh.ry.*( mesh.cell.Dr * G(:,:,i) ) ...
            - mesh.sy.*( mesh.cell.Ds * G(:,:,i) ) ...
            + ( mesh.cell.LIFT * ( mesh.Js .* dFlux(:,:,i) ))./ mesh.J;
    end
end

% for n = 1:obj.Nedge
%     edge = obj.edgeUnion(n);
%     m1 = edge.FToM(1);
%     m2 = edge.FToM(2);
%     mesh1 = obj.meshUnion(m1);
%     mesh2 = obj.meshUnion(m2);
%     for m = 1:edge.M
%         k1 = edge.FToE(1, m);
%         k2 = edge.FToE(2, m);
% 
%         f1 = fphys{m1}(:, k1, :);
%         f2 = fphys{m2}(:, k2, :);
% 
%         % calculate the NToQ1 * f1 and NToQ2 * f2 in each page
%         fq1 = sum( bsxfun( @times, edge.NToQ1(:,:,m), permute(f1, [2,1,3]) ), 2 );
%         fq2 = sum( bsxfun( @times, edge.NToQ2(:,:,m), permute(f2, [2,1,3]) ), 2 );
% 
%         fluxS = obj.matEvaluateEdgeFlux( fq1, fq2, edge.nxq(:,m), edge.nyq(:,m) );
%         [Eq1, Gq1] = obj.matEvaluateFlux( mesh1, fq1 );
%         [Eq2, Gq2] = obj.matEvaluateFlux( mesh2, fq2 );
% 
%         for i = 1:obj.Nvar
%             surfInt1 = edge.NToQ1(:,:,m)' *( edge.Jq(:, m) .* edge.bcell.wq .* ...
%                 ( Eq1(:,:,i).*edge.nxq(:,m) + Gq1.*edge.nyq(:,m) - fluxS(:,:,i) ) );
%             surfInt2 = edge.NToQ2(:,:,m)' *( edge.Jq(:, m) .* edge.bcell.wq .* ...
%                 ( - Eq2(:,:,i).*edge.nxq(:,m) - Gq2.*edge.nyq(:,m) + fluxS(:,:,i) ) );
%             obj.frhs{m1}(:, k1, i) = obj.frhs{m1}(:, k1, i) + ( mesh1.cell.invM * surfInt1 )./mesh1.J(:, k1);
%             obj.frhs{m2}(:, k2, i) = obj.frhs{m2}(:, k2, i) - ( mesh2.cell.invM * surfInt2 )./mesh2.J(:, k2);
%         end
%     end
% end

obj.matEvaluateSourceTerm( fphys );
end