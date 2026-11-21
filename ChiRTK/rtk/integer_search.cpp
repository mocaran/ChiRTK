
#include "gnss_core_types.h"
#include <math.h>
#include <memory.h>

#define kMaxSearchNodes 100000           /* maximum count of search loop */

#define kSign(x)    ((x)<=0.0?-1.0:1.0)
#define kRound(x)   (floor((x)+0.5))
#define kSwap(x,y)  do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)


/*!
@brief 对协方差矩阵做LD分解，得到单位下三角矩阵与对角矩阵
@param n 矩阵阶数，即浮点参数个数
@param Q 输入的对称正定协方差矩阵，n×n，按列主序存储
@param L 输出的单位下三角矩阵，n×n，按列主序存储
@param D 输出的对角矩阵元素，长度n
@return 0=分解成功，-1=对角元素非正说明矩阵非正定
@note 分解式为Q=L'*diag(D)*L；内部复制Q到临时数组，不修改输入矩阵
*/
static int CholeskyLdDecomposition(int n, const double* Q, double* L, double* D)
{
    int i, j, k, info = 0;
    double a, * A;

    A = new double[n * n];
    memcpy(A, Q, sizeof(double) * n * n);
    for (i = n - 1; i >= 0; i--) {
        if ((D[i] = A[i + i * n]) <= 0.0) { info = -1; break; }
        a = sqrt(D[i]);
        for (j = 0; j <= i; j++) L[i + j * n] = A[i + j * n] / a;
        for (j = 0; j <= i - 1; j++) for (k = 0; k <= j; k++) A[j + k * n] -= L[i + k * n] * L[i + j * n];
        for (j = 0; j <= i; j++) L[i + j * n] /= L[i + i * n];
    }

    delete[]A;
    return info;
}


/*!
@brief 对L矩阵与Z矩阵同时做整数高斯变换，消去指定非对角元素
@param n 矩阵阶数
@param L 单位下三角矩阵，n×n，按列主序就地修改
@param Z 整数变换矩阵，n×n，按列主序就地修改
@param i 被消去的行下标
@param j 目标列下标
@return 无
@note 变换系数mu由L[i+j*n]四舍五入得到，mu为0时不修改任何元素
*/
void GaussIntegerTransform(int n, double* L, double* Z, int i, int j)
{
    int k, mu;

    if ((mu = (int)kRound(L[i + j * n])) != 0) {
        for (k = i; k < n; k++) L[k + n * j] -= (double)mu * L[k + i * n];
        for (k = 0; k < n; k++) Z[k + n * j] -= (double)mu * Z[k + i * n];
    }
}


/*!
@brief 交换相邻两列并同步更新L、D与Z，完成一次约化置换
@param n 矩阵阶数
@param L 单位下三角矩阵，n×n，按列主序就地修改
@param D 对角矩阵元素，长度n，就地修改
@param j 参与交换的列下标，与j+1列交换
@param del 由D[j]、D[j+1]与L[j+1+j*n]算出的新对角元素
@param Z 整数变换矩阵，n×n，按列主序就地修改
@return 无
@note 仅在约化满足交换判据时调用，保持Qz分解关系不变
*/
void SwapBasisColumns(int n, double* L, double* D, int j, double del, double* Z)
{
    int k;
    double eta, lam, a0, a1;

    eta = D[j] / del;
    lam = D[j + 1] * L[j + 1 + j * n] / del;
    D[j] = eta * D[j + 1]; D[j + 1] = del;
    for (k = 0; k <= j - 1; k++) {
        a0 = L[j + k * n]; a1 = L[j + 1 + k * n];
        L[j + k * n] = -L[j + 1 + j * n] * a0 + a1;
        L[j + 1 + k * n] = eta * a0 + lam * a1;
    }
    L[j + 1 + j * n] = lam;
    for (k = j + 2; k < n; k++) kSwap(L[k + j * n], L[k + (j + 1) * n]);
    for (k = 0; k < n; k++) kSwap(Z[k + j * n], Z[k + (j + 1) * n]);
}

/*!
@brief 对L、D执行LAMBDA约化，并累积整数变换矩阵Z
@param n 矩阵阶数
@param L 单位下三角矩阵，n×n，按列主序就地修改
@param D 对角矩阵元素，长度n，就地修改
@param Z 输出的整数变换矩阵，n×n，按列主序
@return 无
@note 约化后满足Qz=Z'*Q*Z=L'*diag(D)*L；交换判据加1E-6容差以容忍数值误差
*/
void LllBasisReduction(int n, double* L, double* D, double* Z)
{
    int i, j, k;
    double del;

    j = n - 2; k = n - 2;
    while (j >= 0) {
        if (j <= k) for (i = j + 1; i < n; i++) GaussIntegerTransform(n, L, Z, i, j);
        del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
        if (del + 1E-6 < D[j + 1]) { /* compared considering numerical error */
            SwapBasisColumns(n, L, D, j, del, Z);
            k = j; j = n - 2;
        }
        else j--;
    }
}


/*!
@brief 用MLAMBDA方法搜索m组整数最小二乘模糊度候选解
@param n 浮点参数个数，即模糊度维数
@param m 需要保留的固定解个数
@param L 约化后的单位下三角矩阵，n×n，按列主序
@param D 约化后的对角矩阵元素，长度n
@param zs 约化后的浮点解向量，长度n
@param zn 输出的整数候选解，n×m，按列主序，每列一组解
@param s 输出的候选解残差平方和，长度m，按升序排列
@return 0=搜索成功，-1=搜索循环次数超过kMaxSearchNodes上限
@note 采用深度优先椭球搜索并动态收缩搜索半径，结果最后按残差平方和排序
*/
int SearchIntegerCandidates(int n, int m, const double* L, const double* D,
    const double* zs, double* zn, double* s)
{
    int i, j, k, c, nn = 0, imax = 0;
    double newdist, maxdist = 1E99, y;
    double* S, * dist, * zb, * z, * step;

    S = new double[n * n];
    dist = new double[n];
    zb = new double[n];
    z = new double[n];
    step = new double[n];
    memset(S, 0, n * n * sizeof(double));

    k = n - 1; dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = kRound(zb[k]); y = zb[k] - z[k]; step[k] = kSign(y);
    for (c = 0; c < kMaxSearchNodes; c++)
    {
        newdist = dist[k] + y * y / D[k];
        if (newdist < maxdist)
        {
            if (k != 0)
            {
                dist[--k] = newdist;
                for (i = 0; i <= k; i++) S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
                zb[k] = zs[k] + S[k + k * n];
                z[k] = kRound(zb[k]); y = zb[k] - z[k]; step[k] = kSign(y);
            }
            else {
                if (nn < m)
                {
                    if (nn == 0 || newdist > s[imax]) imax = nn;
                    for (i = 0; i < n; i++) zn[i + nn * n] = z[i];
                    s[nn++] = newdist;
                }
                else {
                    if (newdist < s[imax])
                    {
                        for (i = 0; i < n; i++) zn[i + imax * n] = z[i];
                        s[imax] = newdist;
                        for (i = imax = 0; i < m; i++) if (s[imax] < s[i]) imax = i;
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0]; y = zb[0] - z[0]; step[0] = -step[0] - kSign(step[0]);
            }
        }
        else {
            if (k == n - 1) break;
            else {
                k++;
                z[k] += step[k]; y = zb[k] - z[k]; step[k] = -step[k] - kSign(step[k]);
            }
        }
    }
    for (i = 0; i < m - 1; i++) { /* sort by s */
        for (j = i + 1; j < m; j++) {
            if (s[i] < s[j]) continue;
            kSwap(s[i], s[j]);
            for (k = 0; k < n; k++) kSwap(zn[k + i * n], zn[k + j * n]);
        }
    }

    if (c >= kMaxSearchNodes) {
        delete[]S; delete[]dist; delete[]zb; delete[]z; delete[]step;
        return -1;
    }

    delete[]S; delete[]dist; delete[]zb; delete[]z; delete[]step;
    return 0;
}

/*!
@brief 用LAMBDA方法完成整周模糊度固定，输出固定解与残差
@param n 浮点参数个数，即模糊度维数
@param m 需要固定的解个数
@param a 输入的浮点解参数向量，长度n
@param Q 输入的浮点解协方差矩阵，n×n，按列主序
@param F 输出的固定解，n×m，按列主序，每列一组解
@param s 输出的各组固定解残差平方和，长度m
@return 0=固定成功，-1=维数非法或LD分解、搜索失败
@note 流程为LD分解、LAMBDA约化、MLAMBDA搜索，再用Z的逆变换回原参数空间
*/
int SearchIntegerLeastSquares(int n, int m, const double* a, const double* Q, double* F, double* s)
{
    int i, info;
    double* L, * D, * Z, * ZT, * z, * E;

    if (n <= 0 || m <= 0) return -1;
    L = new double[n * n];
    D = new double[n];
    Z = new double[n * n];
    ZT = new double[n * n];
    z = new double[n];
    E = new double[n * 2];

    memset(L, 0, n * n * sizeof(double));
    memset(Z, 0, n * n * sizeof(double));
    for (i = 0; i < n; i++)  Z[i * (n + 1)] = 1.0;

    /* LD factorization */
    if (!(info = CholeskyLdDecomposition(n, Q, L, D)))
    {
        /* lambda reduction */
        LllBasisReduction(n, L, D, Z);
        MultiplyMatrices(n, n, n, 1, Z, a, z);

        /* mlambda search */
        if (!(info = SearchIntegerCandidates(n, m, L, D, z, E, s)))
        {
            InvertMatrix(n, Z, ZT);
            MultiplyMatrices(n, n, n, 1, ZT, E, F);
            MultiplyMatrices(n, n, n, 1, ZT, E + n, F + n);
        }
    }

    delete[] L; delete[] D; delete[]Z; delete[]ZT; delete[]z; delete[]E;
    return info;
}
