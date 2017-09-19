function [ i,j,k ] = trans_ind( N,id )
%TRANS_IND transform the id-th orthgonal basis function into (i, j, k)

Nquad = (N+1)^2;
k = floor( (id-1)/Nquad ) + 1;
id = mod(id, Nquad);

sk = 1;
for i = 0:N
    for j = 0:N
        if (abs(sk-id)<10e-4)
            return;
        end
        sk = sk+1;
    end
end% for
end

