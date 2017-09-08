# NDG-FEM

NDG-FEM (Nodal Discontinuous Galerkin Finite Element Methods) 是采用Matlab语言建立的节点间断伽辽金有限元函数库，其采用 NDG 方法对微分方程进行离散，具有高阶精度、可利用非结构化网格计算等特点。此外，由于格式具有紧致性，其优点同时包括易于并行化，适用hp自适应加密等优势。

Matlab作为广泛应用于科学计算商业软件，其不仅包含高效矩阵计算函数库，同时还有丰富的可视化工具。根据Matlab语言特性，模型采用无积分方法将方程化为矩阵-向量计算形式，在计算一些无法向量化过程时则采用 mex 混合编程方法加速，并调用 OpenMP 并行模式进一步提高计算效率。

## Features

NDG-FEM模型特点包括

* 采用Matlab面向对象语言编写，将单元、网格、计算物理场各个模块相互分离，并提高代码重用性；
* 模型包括多种标准单元，如线单元（一维），三角形（二维）、四边形单元（二维）；
* 模型包含多种水动力问题求解器，包括对流扩散方程（一维、二维），浅水方程（一维、二维）
* 模型采用矩阵形式计算体积分与面积分，有效提高了计算效率；

## Usages

建立一个简单的对流方程求解器包括以下几个主要步骤；

1. 建立标准单元对象；
调用标准三角形类`ndg_lib/std_cell/tri`，生成二阶三角形对象
```
>> tri = ndg_lib.std_cell.tri(2)

tri =

  tri with properties:

           N: 2
        type: Tri
          Nv: 3
         vol: 2
          vr: [3x1 double]
          vs: [3x1 double]
          vt: [3x1 double]
         Nfv: [2 2 2]
        FToV: [2x3 double]
       Nface: 3
    faceType: [Line    Line    Line]
          Np: 6
           r: [6x1 double]
           s: [6x1 double]
           t: [6x1 double]
           V: [6x6 double]
           M: [6x6 double]
          Dr: [6x6 double]
          Ds: [6x6 double]
          Dt: [6x6 double]
       Fmask: [3x3 double]
         Nfp: [3 3 3]
    Nfptotal: 9
        LIFT: [6x9 double]
```
2. 根据网格数据建立网格对象
网格文件位置`Conv2D/mesh/`，包含文件包括
```
.
├── triangle.edge   # 边界文件
├── triangle.ele    # 单元文件
└── triangle.node   # 顶点文件
```
调用三角形网格类 `ndg_lib/mesh/`生成对应网格对象
```
>> casename = 'Conv2D/mesh/triangle';
>> mesh = ndg_lib.mesh.tri_mesh(tri, casename)

In function read_from_file

finish read node file:
-- Conv2D/mesh/triangle.node
Nv =  556

finish read element file:
-- Conv2D/mesh/triangle.ele
K =  1030

finish read boundary file:
-- Conv2D/mesh/triangle.edge


mesh =

  tri_mesh with properties:

        cell: [1x1 ndg_lib.std_cell.tri]
           K: 1030
          Nv: 556
        EToV: [3x1030 double]
        EToR: [1x1030 double]
       EToBS: [3x1030 uint8]
          vx: [556x1 double]
          vy: [556x1 double]
          vz: [556x1 double]
        EToE: [3x1030 double]
        EToF: [3x1030 double]
           x: [6x1030 double]
           y: [6x1030 double]
           z: [6x1030 double]
          rx: [6x1030 double]
          ry: [6x1030 double]
          rz: [6x1030 double]
          sx: [6x1030 double]
          sy: [6x1030 double]
          sz: [6x1030 double]
          tx: [6x1030 double]
          ty: [6x1030 double]
          tz: [6x1030 double]
           J: [6x1030 double]
        eidM: [9x1030 double]
        eidP: [9x1030 double]
     eidtype: [9x1030 uint8]
    eidfscal: [9x1030 double]
          nx: [9x1030 double]
          ny: [9x1030 double]
          nz: [9x1030 double]
          Js: [9x1030 double]
       Nedge: 1585
       Nnode: 4755
          kM: [1585x1 double]
          kP: [1585x1 double]
          fM: [1585x1 double]
          fP: [1585x1 double]
       ftype: [1585x1 uint8]
         idM: [4755x1 double]
         idP: [4755x1 double]
         fpM: [4755x1 double]
         fpP: [4755x1 double]
       fscal: [4755x1 double]
        fnxM: [4755x1 double]
        fnyM: [4755x1 double]
        fnzM: [4755x1 double]
```
3. 利用三角形网格对象生成二维对流求解器对象
```
>> solver = conv2d_rotation(mesh)

solver =

  conv2d_rotation with properties:

      Nfield: 1
         cfl: 0.3000
       ftime: 2.4000
          dt: 0.0042
           u: [6x1030 double]
           v: [6x1030 double]
        mesh: [1x1 ndg_lib.mesh.tri_mesh]
      f_extQ: [6x1030 double]
    obc_file: []
    out_file: []
         f_Q: [6x1030 double]
```
生成求解器后，可以通过对象的`set_obc_file`与`set_out_file`等方法设定开边界文件与输出文件（NetCDF格式），最终调用`RK45_solve`对方程进行求解。
```
>> tic; solver.RK45_solve; toc
Elapsed time is 3.208880 seconds.
```
调用`draw(1)`查看计算场分布图像
![](http://ww1.sinaimg.cn/large/7a1c18a8ly1ffhvpuqo39j21st1cm7cp.jpg)

## Acknowledgment

This software is inspired by _Nodal Discontinuous Galerkin Methods Algorithms, Analysis, and Applications Texts in Applied Mathematics_ (Hesthaven and Tim Warburton, 2008), please refer to [NodalDG](http://www.caam.rice.edu/~timwar/Book/NodalDG.html) for more details.
