function matEvaluateCellInnerFlux( obj, fphys )

for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    n1 = mesh.cell.n1;
    n2 = mesh.cell.n2;
    
    f1 = fphys{m}(n1, :);
    f2 = fphys{m}(n2, :);
    
    [ flux ] = evaluateNumericalFlux( f1, f2, mesh.locNx, mesh.locNy, obj.u0, obj.v0 );
    
    obj.frhs{m} = mesh.cell.VLIFT * ( flux .* mesh.locEdgeJs );
end

end