function [ mesh ] = uniform_mesh( N, Mx, type )
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
xlim = [0, 75]; ylim = [-15, 15];
wall_bc = ndg_lib.bc_type.SlipWall;
face_type = [wall_bc, wall_bc, wall_bc, wall_bc];

My = ceil( Mx*.4 );
switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.get_std_cell(N, type);
        mesh = ndg_lib.mesh.tri_mesh(cell, 'uniform', ...
            {xlim, ylim, Mx, My, face_type});
        
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.get_std_cell(N, type);
        mesh = ndg_lib.mesh.quad_mesh(cell, 'uniform', ...
            {xlim, ylim, Mx, My, face_type});
end
end

