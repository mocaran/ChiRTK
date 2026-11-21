#include "epoch_processing.h"

#include <cstdio>
#include <cstring>



/*!
@brief 处理一个卡尔曼滤波历元：模糊度状态初始化、时间预测、测量更新与结果输出
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param ekf 卡尔曼滤波状态（跨历元保留）
@param out 输出文件与表头状态
@param stats 会话统计量
@return 无
@note 首次进入或观测间隔大于 5 秒时先退回单点定位并重新初始化滤波；预测或更新失败也会把 IsInit 置 false 并退回单点定位
*/
void ProcessEkfEpoch(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult, KalmanState& ekf, SessionOutputs& out, SessionStats& stats)
{
    // 1、备份当前历元的双差数据
    ekf.CurDDObs = rawdata.DDObs;

    // 观测间隔过长时，旧的模糊度状态不可用，需要重新初始化
    if (ekf.IsInit) {
        const double gap = (rawdata.RovEpk_cur.Time.Week - ekf.Time.Week) * 604800.0
            + rawdata.RovEpk_cur.Time.SecOfWeek - ekf.Time.SecOfWeek;
        if (gap <= 0.0 || gap > 5.0) ekf.IsInit = false;
    }

    // 2、初始化滤波：本历元先输出单点定位结果
    if (!ekf.IsInit) {
        PrintSppFallback(rawdata.RovEpk_cur.Time, baseResult, roverResult, "");

        ekf.CurDDObs = rawdata.DDObs;   // 同步当前双差
        ResetKalmanFilter(&rawdata, &baseResult, &roverResult, &ekf);
        ekf.SDObs = rawdata.SdObs;      // 备份单差顺序
        ekf.DDObs = rawdata.DDObs;      // 备份参考星信息
        return;
    }

    // 3、时间预测；失败则退回单点定位
    if (!PredictKalmanState(&rawdata, &ekf)) {
        ekf.IsInit = false;
        PrintSppFallback(rawdata.RovEpk_cur.Time, baseResult, roverResult, "");
        return;
    }

    // 4、测量更新；失败同样退回单点定位
    if (!CorrectKalmanState(&rawdata, &ekf)) {
        ekf.IsInit = false;
        PrintSppFallback(rawdata.RovEpk_cur.Time, baseResult, roverResult, "()");
        return;
    }

    // 5、把当前单差/双差数据交给下一历元使用
    memcpy(&ekf.SDObs, &rawdata.SdObs, sizeof(ekf.SDObs));
    memcpy(&ekf.DDObs, &rawdata.DDObs, sizeof(ekf.DDObs));

    WriteSppResult(rawdata, roverResult, out.spp, out.sppHeader);

    // 6、计算流动站的 BLH 与相对基站的 ENU 偏差
    double basePos[3] = { 0.0, 0.0, 0.0 };
    double roverFinalBlh[3] = { 0.0, 0.0, 0.0 };
    double dENU[3] = { 0.0, 0.0, 0.0 };

    GeodeticToEcef(rawdata.BaseBest.Pos, basePos, R_WGS84, F_WGS84);   // 基站的 ECEF 坐标
    EcefToGeodetic(ekf.X0, roverFinalBlh, R_WGS84, F_WGS84);           // 流动站的 BLH 坐标

    GeodeticCoord baseBlh;
    baseBlh.longitude = rawdata.BaseBest.Pos[0];
    baseBlh.latitude = rawdata.BaseBest.Pos[1];
    baseBlh.height = rawdata.BaseBest.Pos[2];
    EcefOffsetToEnu(basePos, ekf.X0, &baseBlh, dENU);

    if (!out.ekfHeader) { fprintf(out.ekf, "#Time X Y Z(m) dE dN dU(m) Ratio GPS BDS\n"); out.ekfHeader = true; }
    WriteKalmanRecord(&rawdata, &ekf, dENU, out.ekf);

    // 7、实时模式的控制台输出
    double dxyz[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++) dxyz[i] = ekf.X0[i] - baseResult.Position[i];

    if (ANinfo.IsFileData == 0 && ekf.IsFixed && ekf.IsInit) {
        printf("[固定解] 第 %u 周 %.1f 秒   dX=%+0.8f   dY=%+0.8f   dZ=%+0.8f   GPS=%d   BDS=%d   Ratio=%.3f\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek,
            dxyz[0], dxyz[1], dxyz[2],
            ekf.CurDDObs.DDSatNum[0], ekf.CurDDObs.DDSatNum[1], ekf.ratio);
    }
    if (ANinfo.IsFileData == 0 && !ekf.IsFixed && ekf.IsInit) {
        printf("[浮点解] 第 %u 周 %.1f 秒   dX=%+0.8f   dY=%+0.8f   dZ=%+0.8f   GPS=%d   BDS=%d   Ratio=%.3f\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek,
            dxyz[0], dxyz[1], dxyz[2],
            ekf.CurDDObs.DDSatNum[0], ekf.CurDDObs.DDSatNum[1], ekf.ratio);
    }

    WritePlotRecord(&rawdata, roverFinalBlh, &ekf, roverFinalBlh, out.ekfPlot, 1);

    ++stats.savedEpochs;
    if (ekf.IsFixed) ++stats.fixedEpochs;
    else ++stats.floatEpochs;
}
