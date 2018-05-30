# NDG-FEM

NDG-FEM (Nodal Discontinuous Galerkin Finite Element Methods) 是采用Matlab语言建立的节点间断有限元函数库，其采用 NDG 方法对微分方程进行离散，具有高阶精度、可利用非结构化网格计算等特点。此外，由于格式具有紧致性，其优点同时包括易于并行化，适用hp自适应加密等优势。

Matlab作为广泛应用于科学计算商业软件，其不仅包含高效矩阵计算函数库，同时还有丰富的可视化工具。根据Matlab语言特性，模型采用无积分方法将方程化为矩阵-向量计算形式，在计算一些无法向量化过程时则采用 mex 混合编程方法加速，并调用 OpenMP 并行模式进一步提高计算效率。

## Features

NDG-FEM模型特点包括

* 采用Matlab面向对象语言编写，将单元、网格、计算物理场各个模块相互分离，并提高代码重用性；
* 模型包括多种标准单元，如线单元（一维），三角形（二维）、四边形单元（二维）；
* 模型包含多种水动力问题求解器，包括对流扩散方程（一维、二维），浅水方程（一维、二维）
* 模型采用矩阵形式计算体积分与面积分，有效提高了计算效率；

## Contents

1. NdgCell
基本单元模块。包括线单元`StdPoint`、`StdLine`、`StdTri`、`StdQuad`等。函数`getStdCell`生成给定类型的标准单元对象，而枚举类型`StdCellType`则定义了所有单元对饮的编号。
2. NdgMesh
网格模块。定义了

## Usages

## Acknowledgment

This software is inspired by _Nodal Discontinuous Galerkin Methods Algorithms, Analysis, and Applications Texts in Applied Mathematics_ (Hesthaven and Tim Warburton, 2008), please refer to [NodalDG](http://www.caam.rice.edu/~timwar/Book/NodalDG.html) for more details.
