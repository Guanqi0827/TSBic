
##
.R文件:

整理数据集:

1. Cell clustering:单细胞聚类分析

2. Data preprocessing:数据预处理部分，整理出基因表达量数据矩阵

3. Express velocity:整理出基因表达速率矩阵

4. ConsLR: 3.3中分段多项式函数拟合树型结构数据（全部数据）

5. ConsLR fate: 分段多项式函数拟合树型结构数据（小规模数据集实验的数据）

6. Func corr: 函数型数据相关系数的计算

7. test: 3.3中数据差分前后的假设检验



plot:

8. Plaidmodel: 画小规模数据集经典模型的双聚类结果

9. CatCellFate: 画细胞命运饼图、空心环状图

10. blot scatter: 画基因在表达的分支细胞中表达量的时间序列散点图

11. heatmap: 画部分基因表达 0-1 数据矩阵热图，以及小规模数据集双聚类结果热图

12. GeneEnrichment_all: 基因富集分析



双聚类求解算法:

13. GAB: 小规模数据集寻找第一个双聚类的求解算法

14. FindAllBiclusterings: 小规模数据集寻找全部双聚类的求解算法

15. func_all: 全部数据集找双聚类的某些函数和小规模数据集不一样，需要变动

16. GAB_all: 全部数据集寻找第一个双聚类的求解算法

17. FindAllBiclusterings_all: 全部数据集寻找全部双聚类的求解算法

18. FastFAB: 全部数据集寻找全部双聚类的求解算法（多核加速版）

19. Stop rule: 双聚类个数确定方法
