function f = line_orthgonal_func( ind, r )
f = Polylib.JacobiP(r, 0, 0, ind-1);
end

