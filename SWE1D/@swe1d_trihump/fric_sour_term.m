function [ sf ] = fric_sour_term( obj, f_Q )
%FRIC_SOUR_TERM 一维浅水方程摩阻源项
%   湿单元摩阻项采用 manning 公式计算（湿单元所有节点水深非0），
%   干单元摩阻源项为 0。

sf(:,:,1) = zeros(obj.mesh.cell.Np, obj.mesh.K);
sf(:,:,2) = -obj.gra.*obj.n.*f_Q(:,:,2).*abs(f_Q(:,:,2))./f_Q(:,:,1).^(7/3);

sf(:, ~obj.wetflag, 2) = 0;
end

