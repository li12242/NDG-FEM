# ndg_lib 网格介绍

计算网格对象（各类 mesh 对象）主要包含计算域内各个单元几何信息和单元连接信息，此外还有单元面信息单独存储在 edge 对象内。

mesh 对象基本属性包括

| 属性名 | 说明 | 大小 | 类型 |
| :---: | :---: | :---: | :---: |
| cell  | 标准单元对象 | 1 | `ndg_lib.std_cell` |
| K | 单元个数 | 1 | double |
| Nv | 顶点个数 | 1 | double |
| EToV | 单元内顶点编号 | (cell.Nv,K) | double |
| EToR | 单元类型编号 | K | `ndg_lib.mesh_type` |
| EToBS | 单元面类型 | (cell.Nface,K) | `ndg_lib.bc_type` |
| vx,vy,vz | 顶点序号 | Nv | double |
| J | 节点处雅克比行列式 | (cell.Np,K) | double |
