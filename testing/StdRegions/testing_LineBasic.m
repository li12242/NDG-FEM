function testing_LineBasic
if checkOrthogonal
    error('the Orthogonal of base polynimal is incorrect')
else
    fprintf('checking Orthogonal: success\n\n')
end% if

end% func

function flag = checkOrthogonal
% 判断基函数是否正交
% flag = 1  error
% flag = 0  success

% ref:
% [1] Karniadakis G, Sherwin S. Spectral/hp element methods for computational
%     fluid dynamics[M]. Oxford University Press, 2013. 80-80

nOrder = 4; % 正交基最高次数 p
% [^1] Gauss 积分节点个数 n>= (p^2+3)/2
np = floor((nOrder^2 + 3)/2);
[z,w] = Polylib.zwglj(np);
% [z,w] = JacobiGQ(0, 0, np);
V = zeros(numel(z), nOrder+1);
for j=1:nOrder+1
    V(:,j) = Polylib.JacobiP(z(:), 0, 0, j-1);
end
% V = line.VandMatrix;
s = zeros(nOrder+1, nOrder+1);
for i = 1:nOrder+1
    for j = 1:nOrder+1
        s(i,j) = sum((w(:).*V(:,i)).*V(:,j));
    end% for
end% for
% check result, Matrix s shoule be I
err = (abs(s - eye(nOrder+1)));
flag = any(sum(err)>10e-6);
end% func