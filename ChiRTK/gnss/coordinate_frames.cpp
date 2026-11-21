#include "gnss_core_types.h"


/*!
@brief 将大地坐标（经度、纬度、高程）转换为空间直角坐标
@param geodetic 大地坐标[经度,纬度,高程]（度、度、米）
@param ecef 输出的空间直角坐标[3]（米）
@param R 参考椭球长半轴（米）
@param F 参考椭球扁率
@return 无
@note 先按椭球参数计算卯酉圈半径，再投影到三维直角坐标
*/
void GeodeticToEcef(const double geodetic[3], double ecef[3], const double R, const double F)
{
    const double a = R;                              // 长半轴
    const double f = F;                              // 扁率
    const double b = a * (1.0 - f);                  // 短半轴
    const double e2 = (a * a - b * b) / (a * a);     // 第一偏心率平方

    const double lon = geodetic[0] * Rad;            // 经度 (弧度)
    const double lat = geodetic[1] * Rad;            // 纬度 (弧度)
    const double hgt = geodetic[2];                  // 高程 (米)

    const double sinLat = sin(lat);
    const double cosLat = cos(lat);

    // 卯酉圈曲率半径：N = a / sqrt(1 - e^2 sin^2 B)
    const double N = a / sqrt(1.0 - e2 * sinLat * sinLat);

    // 沿法线方向抬升 hgt 后的投影：(N+h)cosB 为赤道面分量
    const double ring = (N + hgt) * cosLat;

    ecef[0] = ring * cos(lon);
    ecef[1] = ring * sin(lon);
    ecef[2] = (N * (1.0 - e2) + hgt) * sinLat;
}


/*!
@brief 将空间直角坐标转换为大地坐标（纬度、经度、高程）
@param ecef 输入的空间直角坐标[3]（米）
@param geodetic 输出的大地坐标[纬度,经度,高程]（度、度、米），顺序与 GeodeticCoord 的 blh 成员一致
@param R 参考椭球长半轴（米）
@param F 参考椭球扁率
@return 无
@note 纬度采用迭代法求解，收敛阈值1e-12弧度，最多迭代1000次
*/
void EcefToGeodetic(const double ecef[3], double geodetic[3], const double R, const double F)
{
    const double a = R;
    const double b = a * (1 - F);
    const double e2 = (a * a - b * b) / (a * a);

    const double x = ecef[0];
    const double y = ecef[1];
    const double z = ecef[2];
    const double p = sqrt(x * x + y * y);

    // 迭代变量：纬度 B（弧度）与高程 H
    int MaxCount = 1000;
    double Threshold = 1e-12;
    int count = 0;
    double delta = 0.0;
    double H = 0.0;

    double B = atan2(z, p * (1 - e2)); // 初始纬度 (弧度)

    do {
        double sinB = sin(B);
        double N = a / sqrt(1.0 - e2 * sinB * sinB);
        H = p / cos(B) - N;

        double B_new = atan2(z + e2 * (a / sqrt(1.0 - e2 * sinB * sinB)) * sinB, p);
        delta = B - B_new;
        B = B_new;
        count++;
    } while (count < MaxCount && std::fabs(delta) > Threshold);

    if (count >= MaxCount) {
        std::cerr << "警告: 未在" << MaxCount
            << "次迭代内收敛（误差：" << std::fabs(delta)
            << " 弧度）" << std::endl;
    }

    geodetic[0] = B * Deg;           // 纬度 (度)
    geodetic[1] = atan2(y, x) * Deg; // 经度 (度)
    geodetic[2] = H;                 // 高程 (米)
}


/*!
@brief 生成由空间直角坐标到站心ENU坐标的旋转矩阵
@param blh 测站大地坐标[纬度,经度,高程]（度、度、米）
@param Rmat 输出的3×3旋转矩阵，按行存放E、N、U分量
@return 无
@note 矩阵与坐标差向量相乘，得到E、N、U方向的增量
*/
void EnuRotationFromGeodetic(const GeodeticCoord* blh, Matrix3d& Rmat)
{
    const double lat = blh->blh[0] * Rad;
    const double lon = blh->blh[1] * Rad;

    const double sinLat = sin(lat);
    const double cosLat = cos(lat);
    const double sinLon = sin(lon);
    const double cosLon = cos(lon);

    Rmat << -sinLon, cosLon, 0.0,
        -sinLat * cosLon, -sinLat * sinLon, cosLat,
        cosLat * cosLon, cosLat * sinLon, sinLat;
}


/*!
@brief 计算卫星与测站位置差在E、N、U方向的投影
@param Xs 卫星（或目标点）坐标（ECEF坐标，米）
@param Xr 测站坐标（ECEF坐标，米）
@param blh 测站大地坐标[纬度,经度,高程]（度、度、米）
@param dENU 输出的E、N、U方向坐标差（米）
@return 无
@note 内部调用EnuRotationFromGeodetic构造旋转矩阵后作坐标变换
*/
void EcefOffsetToEnu(const double Xs[], const double Xr[], const GeodeticCoord* blh, double dENU[])
{
    Matrix3d Rmat;
    EnuRotationFromGeodetic(blh, Rmat);  // E/N/U方向旋转矩阵

    Vector3d delta;
    delta << Xs[0] - Xr[0], Xs[1] - Xr[1], Xs[2] - Xr[2];

    Vector3d d_enu = Rmat * delta;       // R 为正交阵，直接左乘即得站心分量

    dENU[0] = d_enu(0); // E
    dENU[1] = d_enu(1); // N
    dENU[2] = d_enu(2); // U
}


/*!
@brief 根据测站与卫星坐标计算卫星高度角和方位角
@param Xr 测站坐标（ECEF坐标，米）
@param Xs 卫星坐标（ECEF坐标，米）
@param blh 测站大地坐标[纬度,经度,高程]（度、度、米）
@param Elev 输出的卫星高度角（度，0~90）
@param Azim 输出的卫星方位角（度，0~360，从正北顺时针）
@return 无
@note 方位角以正北为0度、向东为正，负值加360度归化
*/
void SatelliteElevationAzimuth(const double Xr[], const double Xs[], const GeodeticCoord* blh, double* Elev, double* Azim)
{
    double dENU[3] = { 0.0, 0.0, 0.0 };
    EcefOffsetToEnu(Xs, Xr, blh, dENU);

    const double east = dENU[0];
    const double north = dENU[1];
    const double up = dENU[2];
    const double horizontal = sqrt(east * east + north * north);

    *Elev = atan2(up, horizontal) * Deg;
    *Azim = atan2(east, north) * Deg;  // E在前、N在后
    if (*Azim < 0.0) *Azim += 360.0;
}
