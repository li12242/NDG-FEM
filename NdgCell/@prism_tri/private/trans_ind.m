function [ i, j, k ] = trans_ind( N, ind )
%TRANS_IND Summary of this function goes here
%   Detailed explanation goes here

Ntri = (N+1)*(N+2)/2;
k = floor( (ind-1)/Ntri ) + 1;
ind = mod(ind, Ntri);

sk = 1;
for i = 0:N
    for j = 0:(N-i)
        if ( abs(sk-ind)<10e-4 )
            return;
        end
        sk = sk+1;
    end
end% for
end

