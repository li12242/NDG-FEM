function [ mesh ] = uniform_mesh( N, M, type )
%UNIFORM_MESH 构造二维均匀三角形或四边形网格。
%   根据给定基函数阶数和标准单元类型，首先构造标准单元对象 cell，随后调用
%   ndg_utility.uniform_mesh 函数库文件生成对应均匀三角形与四边形网格。
%   均匀网格范围根据参数 xmin/xmax 与 ymin/ymax 确定，上下左右四个边界的
%   类型则由 face_type 指定。
% Input
%   N   - 基函数阶数
%   M   - xy 轴上单元格式
%   type - 标准单元类型：三角形/四边形
% Output:
%   mesh - 网格对象
%
xmin = -100; xmax = 100; 
ymin = -100; ymax = 100;
face_type = [ndg_lib.bc_type.ZeroGrad,...
    ndg_lib.bc_type.ZeroGrad, ...
    ndg_lib.bc_type.ZeroGrad, ...
    ndg_lib.bc_type.ZeroGrad];

switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.ndg_cell(N, type);
        [K,EToV,Nv,VX,VY,EToBS,EToR] = ...
            ndg_utility.uniform_mesh.tri_mesh(M, M, ...
            xmin, xmax, ymin, ymax, face_type);
        mesh = ndg_lib.mesh.tri_mesh(cell, Nv, VX, VY, ...
            K, EToV, EToR, EToBS);
        
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.ndg_cell(N, type);
        [K,EToV,Nv,VX,VY,EToBS,EToR] = ...
            ndg_utility.uniform_mesh.quad_mesh(M, M, ...
            xmin, xmax, ymin, ymax, face_type);
        mesh = ndg_lib.mesh.quad_mesh(cell, Nv, VX, VY, ...
            K, EToV, EToR, EToBS);
end

end

