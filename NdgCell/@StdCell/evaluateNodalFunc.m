function [ func ] = evaluateNodalFunc(obj, r, s, t)
func = zeros(numel(r), obj.Np);
for n = 1:obj.Np
    func(:, n) = obj.evaluateOrthogonalFunc(obj.N, n, r, s, t);
end
func = func/obj.V;
end