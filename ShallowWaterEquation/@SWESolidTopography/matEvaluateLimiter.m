function fphys = matEvaluateLimiter( obj, fphys )
for n = 1:obj.Nvar
    fieldId = obj.varFieldIndex(n);
    fphys = obj.limiter.matVertLimit( fphys, fieldId );
end
end% func