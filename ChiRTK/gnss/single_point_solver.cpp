#include"gnss_core_types.h"


/*!
@brief 利用多普勒观测值进行单点测速最小二乘解算
@param Epoch 当前历元观测数据，含多普勒与卫星速度
@param Result 输入定位坐标，输出速度与速度中误差
@return 无
@note 有效卫星少于4颗或法方程奇异时直接返回，不修改结果；
      观测方程按"有效行"紧凑组装，不再为无效卫星保留全零行
*/
void SolveSinglePointVelocity(EpochObservation* Epoch, SinglePointResult* Result)
{
    Result->Time = Epoch->Time;                  // 记录当前历元时间

    double dirX[MAXCHANNUM], dirY[MAXCHANNUM], dirZ[MAXCHANNUM];
    double residual[MAXCHANNUM];
    int rows = 0;

    for (int i = 0; i < Epoch->SatNum; i++)
    {
        if (!Epoch->SatObs[i].Valid || !Epoch->SatPVT[i].Valid) continue;

        // 接收机指向卫星的视线向量与其模长
        const double lineX = Result->Position[0] - Epoch->SatPVT[i].SatPos[0];
        const double lineY = Result->Position[1] - Epoch->SatPVT[i].SatPos[1];
        const double lineZ = Result->Position[2] - Epoch->SatPVT[i].SatPos[2];
        const double range = sqrt(lineX * lineX + lineY * lineY + lineZ * lineZ);

        if (range < 1.0) continue;               // 几何异常，舍弃该卫星

        // 卫星相对测站的径向速度（多普勒对应量）
        const double rangeRate = -(lineX * Epoch->SatPVT[i].SatVel[0]
            + lineY * Epoch->SatPVT[i].SatVel[1]
            + lineZ * Epoch->SatPVT[i].SatVel[2]) / range;

        // 速度观测方程的系数：视线方向余弦 + 钟漂项
        dirX[rows] = lineX / range;
        dirY[rows] = lineY / range;
        dirZ[rows] = lineZ / range;

        residual[rows] = Epoch->SatObs[i].D[0]
            - (rangeRate - C_Light * Epoch->SatPVT[i].SatClkSft);

        rows++;
    }

    if (rows < 4) return;                        // 观测不足，无法解算

    // 紧凑组装观测矩阵与残差向量
    MatrixXd B = MatrixXd::Zero(rows, 4);
    VectorXd w(rows);

    for (int r = 0; r < rows; r++)
    {
        B(r, 0) = dirX[r];
        B(r, 1) = dirY[r];
        B(r, 2) = dirZ[r];
        B(r, 3) = 1.0;
        w(r) = residual[r];
    }

    // 法方程 N·X = W
    const MatrixXd N = B.transpose() * B;
    const VectorXd W = B.transpose() * w;

    // 用全主元 LU 判定可逆性，替代行列式等于零的判断
    FullPivLU<MatrixXd> lu(N);
    if (!lu.isInvertible()) return;

    const VectorXd X = lu.solve(W);
    if (!X.allFinite()) return;

    Result->Velocity[0] = X(0);
    Result->Velocity[1] = X(1);
    Result->Velocity[2] = X(2);
}


/*!
@brief 采用伪距IF组合进行GPS/BDS单点定位最小二乘解算
@param Epoch 当前历元观测数据，解算后更新SatPVT
@param Raw 原始数据，提供GPS与BDS广播星历
@param Result 输出的定位结果，含坐标、钟差、PDOP与中误差
@return true=解算成功，false=有效卫星数不足
@note 观测矩阵按实际参与的系统现场确定列数（三坐标+每系统一个钟差），
      不再先凑成五列再裁剪；迭代收敛阈值 1e-7，最多 16 次
*/
bool SolveSinglePointPosition(EpochObservation* Epoch, SolverData* Raw, SinglePointResult* Result)
{
    Result->Time = Epoch->Time;

    Vector3d pos(Result->Position[0], Result->Position[1], Result->Position[2]);  // 初值取上一历元
    double clkGps = 0.0;                          // GPS 接收机钟差（米）
    double clkBds = 0.0;                          // BDS 接收机钟差（米）

    bool converged = false;
    int iteration = 0;

    while (!converged && iteration < 16)
    {
        // 1、按当前近似坐标重算各卫星的信号发射时刻状态
        SatelliteStateAtSignalTime(Epoch, Raw->GpsEph, Raw->BdsEph, pos.data());

        // 2、逐卫星组装观测方程，只保留有效行
        double coefX[MAXCHANNUM], coefY[MAXCHANNUM], coefZ[MAXCHANNUM];
        double residual[MAXCHANNUM];
        int sysTag[MAXCHANNUM];                   // 0=GPS，1=BDS
        int rows = 0;
        int gpsCount = 0;
        int bdsCount = 0;

        for (int i = 0; i < Epoch->SatNum; i++)
        {
            if (!Epoch->SatObs[i].Valid || !Epoch->SatPVT[i].Valid) continue;

            const bool isGps = (Epoch->SatObs[i].System == GPS);
            const bool isBds = (Epoch->SatObs[i].System == BDS);
            if (!isGps && !isBds) continue;

            // 接收机到卫星的距离
            const double dx = Epoch->SatPVT[i].SatPos[0] - pos[0];
            const double dy = Epoch->SatPVT[i].SatPos[1] - pos[1];
            const double dz = Epoch->SatPVT[i].SatPos[2] - pos[2];
            const double range = sqrt(dx * dx + dy * dy + dz * dz);

            if (range < 1.0) continue;

            // 方向余弦取负号：对应"坐标改正量"的偏导
            coefX[rows] = -dx / range;
            coefY[rows] = -dy / range;
            coefZ[rows] = -dz / range;

            // 模型值 = 几何距离 + 钟差 - 卫星钟差·c + 对流层 + （BDS）群延迟
            double model = range + Epoch->SatPVT[i].TropCorr
                - C_Light * Epoch->SatPVT[i].SatClkOft;

            if (isGps)
            {
                model += clkGps;
                sysTag[rows] = 0;
                gpsCount++;
            }
            else
            {
                model += clkBds;
                model += (FG1_BDS * FG1_BDS * C_Light * Epoch->SatPVT[i].Tgd1)
                    / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);
                sysTag[rows] = 1;
                bdsCount++;
            }

            residual[rows] = Epoch->ComObs[i].PIF - model;
            rows++;
        }

        // 3、参数个数 = 3 个坐标 + 每个参与系统 1 个钟差
        const int paramCount = 3 + (gpsCount > 0 ? 1 : 0) + (bdsCount > 0 ? 1 : 0);
        if (rows < paramCount) return false;

        const int gpsCol = (gpsCount > 0) ? 3 : -1;
        const int bdsCol = (bdsCount > 0) ? (gpsCount > 0 ? 4 : 3) : -1;

        // 4、组装紧凑的观测矩阵与残差向量
        MatrixXd B = MatrixXd::Zero(rows, paramCount);
        VectorXd w(rows);

        for (int r = 0; r < rows; r++)
        {
            B(r, 0) = coefX[r];
            B(r, 1) = coefY[r];
            B(r, 2) = coefZ[r];

            if (sysTag[r] == 0 && gpsCol >= 0) B(r, gpsCol) = 1.0;
            if (sysTag[r] == 1 && bdsCol >= 0) B(r, bdsCol) = 1.0;

            w(r) = residual[r];
        }

        // 5、等权最小二乘：N·δ = Bᵀ·w
        const MatrixXd N = B.transpose() * B;
        const VectorXd rhs = B.transpose() * w;

        LDLT<MatrixXd> ldlt(N);
        if (ldlt.info() != Success) return false;

        const VectorXd delta = ldlt.solve(rhs);
        if (!delta.allFinite()) return false;

        // 6、收敛判据：参数改正量范数
        converged = (delta.norm() < 1e-7);

        // 7、参数更新
        pos += delta.head<3>();
        if (gpsCol >= 0) clkGps += delta[gpsCol];
        if (bdsCol >= 0) clkBds += delta[bdsCol];

        // 8、精度指标：由法方程逆阵的对角元给 PDOP
        const MatrixXd Ninv = N.inverse();
        if (Ninv.allFinite())
        {
            Result->PDOP = sqrt(Ninv(0, 0) + Ninv(1, 1) + Ninv(2, 2));
        }

        // 9、验后单位权中误差
        const VectorXd v = B * delta - w;
        Result->SigmaPos = sqrt(v.dot(v) / (rows - paramCount));

        Result->GPSSatNum = gpsCount;
        Result->BDSSatNum = bdsCount;
        Result->AllSatNum = rows;

        iteration++;
    }

    // 10、保存最终结果
    Result->RcvClkOft[0] = clkGps;
    Result->RcvClkOft[1] = clkBds;

    Result->Position[0] = pos[0];
    Result->Position[1] = pos[1];
    Result->Position[2] = pos[2];

    return true;
}


/*!
@brief 计算各卫星信号发射时刻的位置、速度、钟差及对流层改正
@param Epk 当前历元观测数据，输入观测值并输出SatPVT
@param GPSEph GPS星历数组首地址
@param BDSEph BDS星历数组首地址
@param UserPos 接收机近似位置（ECEF坐标，米）
@return 无
@note 先迭代卫星钟差，再作地球自转改正，并计算高度角、方位角与对流层延迟；
      GPS 与 BDS 的差异集中在地球自转速率与椭球参数上，故只保留一处公共流程
*/
void SatelliteStateAtSignalTime(EpochObservation* Epk, EphemerisRecord* GPSEph, EphemerisRecord* BDSEph, double UserPos[3])
{
    const GpsTime epochTime = Epk->Time;

    for (int i = 0; i < Epk->SatNum; i++)
    {
        SatelliteState* state = Epk->SatPVT + i;

        // 1、先判定星历可获得性，过期则整颗卫星置无效
        if (!SatelliteClockCorrection(Epk->SatObs[i].Prn, Epk->SatObs[i].System,
            &Epk->Time, GPSEph, BDSEph, state))
        {
            state->Valid = false;
            continue;
        }

        const GnssSystem sys = Epk->SatObs[i].System;
        if (sys != GPS && sys != BDS)
        {
            state->Valid = false;
            continue;
        }

        const bool isGps = (sys == GPS);
        EphemerisRecord* eph = isGps ? (GPSEph + Epk->SatObs[i].Prn - 1)
            : (BDSEph + Epk->SatObs[i].Prn - 1);

        // 2、以钟差为未知量迭代求信号发射时刻
        const double pseudorange = Epk->SatObs[i].P[0];
        GpsTime sigTime;
        double clkOffset = 0.0;

        for (int iter = 0; iter < 16; iter++)
        {
            sigTime.Week = epochTime.Week;
            sigTime.SecOfWeek = epochTime.SecOfWeek - pseudorange / C_Light - clkOffset;

            // BDS 时间系统相对 GPS 有 1356 周与 14 秒的偏移
            const double sysWeek = isGps ? sigTime.Week : (sigTime.Week - 1356.0);
            const double sysSec = isGps ? sigTime.SecOfWeek : (sigTime.SecOfWeek - 14.0);

            const double ageOfClock = (sysWeek - eph->TOC.Week) * 604800.0
                + sysSec - eph->TOC.SecOfWeek;

            const double newOffset = eph->ClkBias + eph->ClkDrift * ageOfClock
                + eph->ClkDriftRate * ageOfClock * ageOfClock;

            const double change = fabs(newOffset - clkOffset);
            clkOffset = newOffset;

            if (change <= 1e-12) break;
        }

        state->SatClkOft = clkOffset;

        // 3、按发射时刻计算卫星位置、速度、钟差与钟速
        if (isGps)
        {
            ComputeGpsSatelliteState(Epk->SatObs[i].Prn, &sigTime, eph, state);
        }
        else
        {
            ComputeBdsSatelliteState(Epk->SatObs[i].Prn, &sigTime, eph, state);
        }

        // 4、信号传播时间与地球自转改正
        const double dx = state->SatPos[0] - UserPos[0];
        const double dy = state->SatPos[1] - UserPos[1];
        const double dz = state->SatPos[2] - UserPos[2];
        const double travelTime = sqrt(dx * dx + dy * dy + dz * dz) / C_Light;

        const double rate = isGps ? W_WGS84 : W_CGCS2000;
        const double alpha = rate * travelTime;         // 传播期间地球转过的角度
        const double cosA = cos(alpha);
        const double sinA = sin(alpha);

        // 先把原值取出来，避免就地覆盖后互相污染
        const double px = state->SatPos[0];
        const double py = state->SatPos[1];
        const double pz = state->SatPos[2];
        const double vx = state->SatVel[0];
        const double vy = state->SatVel[1];
        const double vz = state->SatVel[2];

        state->SatPos[0] = cosA * px + sinA * py;        // 绕 z 轴顺时针转过 alpha
        state->SatPos[1] = -sinA * px + cosA * py;
        state->SatPos[2] = pz;                           // z 分量在旋转下不变

        state->SatVel[0] = cosA * vx + sinA * vy;
        state->SatVel[1] = -sinA * vx + cosA * vy;
        state->SatVel[2] = vz;

        // 5、高度角、方位角与对流层延迟
        GeodeticCoord siteGeodetic;
        const double axisRadius = isGps ? R_WGS84 : R_CGCS2000;
        const double flattening = isGps ? F_WGS84 : F_CGCS2000;

        EcefToGeodetic(UserPos, siteGeodetic.blh, axisRadius, flattening);
        SatelliteElevationAzimuth(UserPos, state->SatPos, &siteGeodetic,
            &state->Elevation, &state->Azimuth);

        state->TropCorr = HopfieldTroposphereDelay(siteGeodetic.blh[2], state->Elevation);
    }
}
