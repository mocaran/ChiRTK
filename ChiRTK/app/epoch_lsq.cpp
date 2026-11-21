#include "epoch_processing.h"

#include <cstdio>



/*!
@brief 处理一个最小二乘历元：浮点解、整数固定、固定解与结果输出
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param floatSolution 浮点解结果（跨历元保留）
@param ekf 卡尔曼滤波状态，仅用于绘图格式参数（最小二乘模式不使用其中数值）
@param out 输出文件与表头状态
@param stats 会话统计量
@return 无
@note 浮点解失败时退回单点定位并往绘图文件写一条 SPP 记录；Ratio 检验未通过或固定解失败时本历元不输出
*/
void ProcessLsqEpoch(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult, FloatSolution& floatSolution, KalmanState& ekf,
    SessionOutputs& out, SessionStats& stats)
{
    // 1、浮点解；失败则退回单点定位，并往绘图文件写一条 SPP 记录
    if (!EstimateFloatSolution(&rawdata, &baseResult, &roverResult, &floatSolution)) {
        PrintSppFallback(rawdata.RovEpk_cur.Time, baseResult, roverResult, "");

        double sppBlh[3];
        EcefToGeodetic(roverResult.Position, sppBlh, R_WGS84, F_WGS84);
        fprintf(out.lsqPlot, "%u %.f %0.18f %0.18f %0.18f 0 3\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek,
            sppBlh[0], sppBlh[1], sppBlh[2]);
        ++stats.savedEpochs;
        ++stats.sppEpochs;
        return;
    }

    // 2、整数最小二乘固定；失败时本历元不输出
    if (SearchIntegerLeastSquares(rawdata.DDObs.Sats * 2, 2, floatSolution.N, floatSolution.Qn,
        rawdata.DDObs.FixedAmb, rawdata.DDObs.ResAmb) != 0)
        return;

    // 3、用次优/最优残差平方和之比做 Ratio 检验
    rawdata.DDObs.Ratio = rawdata.DDObs.ResAmb[1] / rawdata.DDObs.ResAmb[0];
    if (rawdata.DDObs.Ratio >= ANinfo.RatioThres) {
        if (rawdata.DDObs.Ratio >= 100) rawdata.DDObs.Ratio = 100;   // 大于 100 时按 100 记
        rawdata.DDObs.bFixed = true;
    }
    else {
        rawdata.DDObs.bFixed = false;
    }

    // 4、固定解基线；失败时本历元不输出
    if (!EstimateFixedBaseline(&rawdata, &baseResult, &roverResult, &floatSolution)) {
        printf("第 %u 周 %.1f 秒：固定解计算失败，本历元不输出结果。\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek);
        return;
    }

    WriteSppResult(rawdata, roverResult, out.spp, out.sppHeader);

    // 5、基站坐标只算一次，浮点解和固定解都用它作 ENU 的参考点
    double basePos[3] = { 0.0, 0.0, 0.0 };
    GeodeticToEcef(rawdata.BaseBest.Pos, basePos, R_WGS84, F_WGS84);

    GeodeticCoord baseBlh;
    baseBlh.longitude = rawdata.BaseBest.Pos[0];
    baseBlh.latitude = rawdata.BaseBest.Pos[1];
    baseBlh.height = rawdata.BaseBest.Pos[2];

    // 6、浮点解输出
    double roverFloatXYZ[3] = { 0.0, 0.0, 0.0 };
    double roverFloatBlh[3] = { 0.0, 0.0, 0.0 };
    double dENUFloat[3] = { 0.0, 0.0, 0.0 };
    for (int j = 0; j < 3; j++) roverFloatXYZ[j] = floatSolution.dX[j] + basePos[j];
    EcefToGeodetic(roverFloatXYZ, roverFloatBlh, R_WGS84, F_WGS84);
    EcefOffsetToEnu(basePos, roverFloatXYZ, &baseBlh, dENUFloat);

    if (!out.floatHeader) {
        fprintf(out.floatSolution, "#Time X Y Z(m) dE dN dU(m) Sigma PDOP HDOP VDOP GPS BDS\n");
        out.floatHeader = true;
    }
    WriteFloatRecord(&rawdata, roverFloatXYZ, dENUFloat, &floatSolution, out.floatSolution);

    // 7、固定解输出
    double roverFixedXYZ[3] = { 0.0, 0.0, 0.0 };
    double roverFixedBlh[3] = { 0.0, 0.0, 0.0 };
    double dENUFixed[3] = { 0.0, 0.0, 0.0 };
    for (int j = 0; j < 3; j++) roverFixedXYZ[j] = rawdata.DDObs.dPos[j] + basePos[j];
    EcefToGeodetic(roverFixedXYZ, roverFixedBlh, R_WGS84, F_WGS84);
    EcefOffsetToEnu(basePos, roverFixedXYZ, &baseBlh, dENUFixed);

    if (!out.fixedHeader) {
        fprintf(out.fixedSolution, "#Week Sec X Y Z(m) dE dN dU(m) Ratio GPS BDS\n");
        out.fixedHeader = true;
    }
    WriteFixedRecord(&rawdata, roverFixedXYZ, dENUFixed, out.fixedSolution);


    // 8、实时模式的控制台输出
    if (ANinfo.IsFileData == 0 && rawdata.DDObs.bFixed) {
        printf("[固定解] 第 %u 周 %.1f 秒   dX=%+0.8f   dY=%+0.8f   dZ=%+0.8f   PDOP=%.3f   GPS=%d   BDS=%d   Ratio=%.3f\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek,
            rawdata.DDObs.dPos[0], rawdata.DDObs.dPos[1], rawdata.DDObs.dPos[2],
            floatSolution.PDOP, rawdata.DDObs.DDSatNum[0], rawdata.DDObs.DDSatNum[1],
            rawdata.DDObs.Ratio);
    }
    if (ANinfo.IsFileData == 0 && !rawdata.DDObs.bFixed) {
        printf("[浮点解] 第 %u 周 %.1f 秒   dX=%+0.8f   dY=%+0.8f   dZ=%+0.8f   PDOP=%.3f   GPS=%d   BDS=%d   Ratio=%.3f\n",
            rawdata.RovEpk_cur.Time.Week, rawdata.RovEpk_cur.Time.SecOfWeek,
            floatSolution.dX[0], floatSolution.dX[1], floatSolution.dX[2],
            floatSolution.PDOP, rawdata.DDObs.DDSatNum[0], rawdata.DDObs.DDSatNum[1],
            rawdata.DDObs.Ratio);
    }

    WritePlotRecord(&rawdata, roverFixedBlh, &ekf, roverFixedBlh, out.lsqPlot, 2);

    ++stats.savedEpochs;
    if (rawdata.DDObs.bFixed) ++stats.fixedEpochs;
    else ++stats.floatEpochs;
}
