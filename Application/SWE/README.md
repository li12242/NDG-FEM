\page SWE Shallow Water Equations

# 1. 浅水方程介绍

浅水方程描述的是关于质量与动量守恒方程，其表达式为

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U}),
$$

其中 $\mathbf{U} = (h, hu, hv)$ 为守恒变量，$\mathbf{F}(\mathbf{U})$ 与 $\mathbf{S}(\mathbf{U})$ 分别为通量项与源项，表达式为

$$
\mathbf{F}(\mathbf{U}, z) = \left[ 
\begin{array}{cc} hu & hv \cr 
hu^2 + \frac{1}{2}g h^2 & huv \cr 
huv & hv^2 + \frac{1}{2}g h^2
\end{array} \right],
$$
$$
\mathbf{S}(\mathbf{U}, z) = \left[ \begin{array}{c} 0 \\ 
-g h \frac{\partial z}{\partial x} + \frac{u n^2 \sqrt{u^2 + v^2}}{h^{4/3}} \cr 
-g h \frac{\partial z}{\partial y} + \frac{v n^2 \sqrt{u^2 + v^2}}{h^{4/3}} \end{array} \right].
$$

当源项 $\mathbf{S}(\mathbf{U}) = 0$ 时，方程为守恒方程形式，可采用标准 NDG 格式进行离散求解。
对于源项，我们需要一些特殊处理方法以保证计算精度和效率的平衡。首先将源项用基函数近似

$$
\mathbf{S}(\mathbf{U}) \approx \mathbf{S}_h = \sum_{j=1}^{N_p} \mathbf{S}_j \varphi_j(x),
$$

其中，源项数值解 $\mathbf{S}_h$ 系数表达式 $\mathbf{S}_j$ 可以用节点处守恒变量解 $\mathbf{U}_j$ 计算得到。
方程随后乘以试验函数并在控制体内进行积分，这时，根据质量矩阵表达式，可以得出源项积分形式为

$$
\int_{\Omega_k} \varphi_i \mathbf{S}_h \mathrm{d}A = \sum_{j=1}^{N_p} \left( \int_{\Omega_k} \varphi_i \varphi_j  \mathrm{d}A \right) \mathbf{S}_j = \sum_{j=1}^{N_p} \mathcal{M}_{i,j} \mathbf{S}_j.
$$

最终方程乘以质量矩阵逆矩阵，得到右端项表达式为 

$$
\mathcal{L}(\mathbf{U}_h) = - \mathcal{D}_x \cdot \mathbf{E} - \mathcal{D}_y \cdot \mathbf{G} + \mathcal{M}_{lift} (\mathbf{F} - \mathbf{F}^*) \cdot \mathbf{n} + \mathbf{S}.
$$

# 2. Well-balancing 与干湿处理

在浅水方程模拟中，数值方法满足 well-balancing 性质是非常重要的，其表示数值格式能够在不同地形下准确模拟静水情形

$$
\eta = 0, \quad \mathbf{q} = 0,
$$

我们采用 Kesserwani & Liang (2010), Duran & Marche (2004) 所采用的 pre-balanced SWE 进行求解。其待求解的守恒变量不变，通量项与底坡源项表达式为

$$
\mathbf{F}(\mathbf{U}, z) = \left[ 
\begin{array}{cc} hu & hv \cr 
hu^2 + \frac{1}{2}g( h^2 - z^2 ) & huv \cr 
huv & hv^2 + \frac{1}{2}g( h^2 - z^2 )
\end{array} \right], 
$$
$$
\mathbf{S}(\mathbf{U}, z) = \left[ \begin{array}{c} 0 \\ 
-g \eta \frac{\partial z}{\partial x} + \frac{u n^2 \sqrt{u^2 + v^2}}{h^{4/3}} \cr 
-g \eta \frac{\partial z}{\partial y} + \frac{v n^2 \sqrt{u^2 + v^2}}{h^{4/3}} \end{array} \right].
$$

在浅水方程中，干湿处理也是非常复杂的一个问题。在采用固定网格实际模拟过程中，在落潮等时刻水体开始流出计算单元，此时会导致计算单元水深降为0以下的负值。因此需要对模型进行数值上特殊处理，去除负水深出现的单元，并且保证全局质量守恒。

首先需要将单元分为干单元，湿单元和半干半湿单元三类。
干与湿单元可以通过单元内所有节点节点水深来确定，为了保证干湿方法鲁棒性，我们采用一个大于零的水深阈值 $h_{min} \ge 0$ 来判断节点水深。
当且仅当单元内存在部分节点水深大于 $h_{min} \ge 0$ 时，单元为半干半湿单元。

对于干单元，其单元内由于水深为0，因此没有静水压力与动量，所以将流量 $(hu, hv)$ 强制赋值为零，且通量项 $\mathbf{F} = (\mathbf{E}, \mathbf{G})$ 和源项 $\mathbf{S}$ 也都为零。

对于半干半湿单元，需要对其水深（流量）进行修正。
第一种方法采用 Xing et al. (2010) 的 positivity-preserving 算子对半干半湿单元进行修正，但是这要求所有单元平均水深为非负值

$$
\bar{h}_k = \int_{\Omega_k} h_h \mathrm{d}A \ge 0,
$$

我们采用一种简单的干湿处理方法：对于平均水深小于零的半干半湿单元，强制其平均水深为零。这不能保证水体质量守恒，但是由于半干半湿单元负水深很小，所以对计算结果影响并不大。

另外一种方法为 $p$ 阶自适应方法，其中半干半湿单元采用线性解近似。对于水深为负单元，首先将水深设为0，然后将增加的水体积从所有湿单元中去除，最终保证质量守恒。

# 3. 求解器模块组成

求解器根据控制方程不同可以分为两类

1. 传统 SWE；
2. Pre-balanced SWE；

针对干湿方法也可分为两类

1. 简单干湿处理方法；
2. $p$ 阶自适应干湿方法；

因此，SWE求解器细分应该有4个。

在采用面向对象形式的 NDG-FEM 中，实现这种求解器特化方法大致可分为两种。

## 3.1. 求解器模块实现方法

第一种是通过继承与重载实现求解器特化。这种方法中我们针对每一种求解器需生成一个对象，如干湿判断，well-balance，数值通量等等。
对于其他小型修改，可以令新对象继承自此对象，并针对各自特性对函数进行重载。  
这种方法的缺点是会生成繁多的 SWE 求解器，在定义每个算例时需要根据算例设置修改继承不同求解器父类，此时，算例子类在计算时到底在调用哪个的方法并不清晰。另外，Matlab 子类继承多个父类时初始化也存在的问题，如多个父类中不能重复定义相同名称的属性，还有对象初始化时调用父类构造函数顺序等问题。

第二种方法是通过属性初始化实现求解器特化。这种方法是在 SWE 求解器中声明一个属性，专门实现不同 SWE 求解器中计算特化部分。  
如 Covectional SWE 和 Pre-balanced SWE 两类求解器中，其主要区别在于通量项和底坡源项的计算。因此我们可以在求解器内声明一个 fluxSolver 属性，其储存一个 SWEFluxSolver 对象，这个对象包含 evaluateAdvFluxTerm 和 evaluateBottomSourceTerm 方法，分别负责通量项和底坡源项计算。但是在通量项中不仅包含 well-balanced 性质修改，还可能包含干湿等内容修正，这就引起。   
这种方法缺点是需要在算例中定义不同的参数，并在初始化时进行各求解项的构造，此外不同属性所包含的方法之间可能会出现重叠。

## 3.2. 求解器介绍

综合考虑，模型将采用第一种方法来实现求解器构造，求解器列表与继承关系如下图所示



# 4. 一维浅水方程验证与应用

# 5. 二维浅水方程验证与应用