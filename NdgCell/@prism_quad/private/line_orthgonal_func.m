function f = line_orthgonal_func( ind, r )
% LINE_ORTHGONAL_FUNC
f = Polylib.JacobiP(r, 0, 0, ind-1);
end

