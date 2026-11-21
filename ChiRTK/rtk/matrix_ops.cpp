#include "gnss_core_types.h"

/*!
@brief 用高斯-约当消元法对n阶方阵求逆
@param n 矩阵阶数
@param a 输入的方阵，n×n，按行优先存储，函数内不修改
@param b 输出的逆矩阵，n×n，按行优先存储
@return 1=求逆成功，0=维数非法或矩阵奇异（主元素小于DBL_EPSILON）
@note 每步选取绝对值最大的主元素并记录行列交换，最后按逆序恢复行列次序
*/
int InvertMatrix(int n, double a[], double b[])
{
    int i, j, k, l, u, v, is[200], js[200];
    double d, p;

    if (n <= 0)
    {
        printf("矩阵求逆失败：矩阵维数不匹配。\n");
        return 0;
    }

    // 将输入矩阵赋值给输出矩阵b，下面对b矩阵求逆，a矩阵不变
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            b[i * n + j] = a[i * n + j];
        }
    }

    for (k = 0; k < n; k++)
    {
        d = 0.0;
        for (i = k; i < n; i++)   // 查找右下角方阵中主元素的位置
        {
            for (j = k; j < n; j++)
            {
                l = n * i + j;
                p = fabs(b[l]);
                if (p > d)
                {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }

        if (d < DBL_EPSILON)   // 主元素接近于0，矩阵不可逆
        {
            printf("矩阵求逆失败：出现除以 0。\n");
            return 0;
        }

        if (is[k] != k)  // 对主元素所在的行与右下角方阵的首行进行调换
        {
            for (j = 0; j < n; j++)
            {
                u = k * n + j;
                v = is[k] * n + j;
                p = b[u];
                b[u] = b[v];
                b[v] = p;
            }
        }

        if (js[k] != k)  // 对主元素所在的列与右下角方阵的首列进行调换
        {
            for (i = 0; i < n; i++)
            {
                u = i * n + k;
                v = i * n + js[k];
                p = b[u];
                b[u] = b[v];
                b[v] = p;
            }
        }

        l = k * n + k;
        b[l] = 1.0 / b[l];  // 初等行变换
        for (j = 0; j < n; j++)
        {
            if (j != k)
            {
                u = k * n + j;
                b[u] = b[u] * b[l];
            }
        }
        for (i = 0; i < n; i++)
        {
            if (i != k)
            {
                for (j = 0; j < n; j++)
                {
                    if (j != k)
                    {
                        u = i * n + j;
                        b[u] = b[u] - b[i * n + k] * b[k * n + j];
                    }
                }
            }
        }
        for (i = 0; i < n; i++)
        {
            if (i != k)
            {
                u = i * n + k;
                b[u] = -b[u] * b[l];
            }
        }
    }

    for (k = n - 1; k >= 0; k--)  // 将上面的行列调换重新恢复
    {
        if (js[k] != k)
        {
            for (j = 0; j < n; j++)
            {
                u = k * n + j;
                v = js[k] * n + j;
                p = b[u];
                b[u] = b[v];
                b[v] = p;
            }
        }
        if (is[k] != k)
        {
            for (i = 0; i < n; i++)
            {
                u = i * n + k;
                v = is[k] + i * n;
                p = b[u];
                b[u] = b[v];
                b[v] = p;
            }
        }
    }

    return (1);
}

/*!
@brief 按行优先完成两个矩阵的乘法M3=M1*M2
@param m1 左矩阵行数
@param n1 左矩阵列数，必须等于m2
@param m2 右矩阵行数
@param n2 右矩阵列数
@param M1 左矩阵，m1×n1，按行优先存储
@param M2 右矩阵，m2×n2，按行优先存储
@param M3 输出的结果矩阵，m1×n2，按行优先存储
@return true=相乘成功，false=维数非法
@note 维数检查要求各维为正且n1==m2；结果元素先清零再逐项累加
*/
bool MultiplyMatrices(int m1, int n1, int m2, int n2, const double M1[], const double M2[], double M3[])
{
    // 首先判断矩阵相乘的大小是否合法
    if (m1 <= 0 || n1 <= 0 || m2 <= 0 || n2 <= 0 || n1 != m2)
    {
        printf("矩阵相乘失败：矩阵维数不匹配。\n");
        return false; // 矩阵大小不合法，无法进行矩阵相乘
    }

    // 执行矩阵相乘运算
    for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            M3[i * n2 + j] = 0; // 初始化结果矩阵元素为0
            for (int k = 0; k < n1; k++)
            {
                M3[i * n2 + j] = M1[i * n1 + k] * M2[k * n2 + j] + M3[i * n2 + j]; // 矩阵相乘
            }
        }
    }

    return true;
}
