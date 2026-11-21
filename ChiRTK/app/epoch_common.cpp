#include "epoch_processing.h"

#include <cstdio>



/*!
@brief 同步一个历元的基站/流动站观测数据
@param baseFile 基站观测文件（文件模式）
@param roverFile 流动站观测文件（文件模式）
@param baseSocket 基站数据流（实时模式）
@param roverSocket 流动站数据流（实时模式）
@param rawdata 输出的原始观测数据
@return 1=同步成功，0=时间不齐需要跳过，-1=文件结束，-2=数据流失败
@note 按 ANinfo.IsFileData 在文件同步与实时同步之间二选一
*/
int SynchronizeEpoch(FILE* baseFile, FILE* roverFile, SOCKET& baseSocket, SOCKET& roverSocket,
    SolverData& rawdata)
{
    if (ANinfo.IsFileData == 1) return SyncEpochFromFiles(baseFile, roverFile, &rawdata);
    if (ANinfo.IsFileData == 0) return SyncEpochFromSockets(baseSocket, roverSocket, &rawdata);
    return 0;   // 既不是文件也不是实时模式，本历元不取数
}

/*!
@brief 分别对基站、流动站做单点定位与测速（卫星位置也在这里算出）
@param rawdata 原始观测数据，函数会更新其中的单点定位状态
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@return 无
@note 定位前先调用 RejectObservationOutliers 剔除粗差；只有单点定位成功才继续测速
*/
void SolveStationPositions(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult)
{
    // 基站：先剔除粗差，再单点定位；定位成功才继续测速
    RejectObservationOutliers(&rawdata.BasEpk_cur);
    if (SolveSinglePointPosition(&rawdata.BasEpk_cur, &rawdata, &baseResult))
        SolveSinglePointVelocity(&rawdata.BasEpk_cur, &baseResult);

    // 流动站：同样的步骤
    RejectObservationOutliers(&rawdata.RovEpk_cur);
    if (SolveSinglePointPosition(&rawdata.RovEpk_cur, &rawdata, &roverResult))
        SolveSinglePointVelocity(&rawdata.RovEpk_cur, &roverResult);
}

/*!
@brief 完成进入 RTK 解算前的公共步骤：有效性标记、站间单差、周跳探测、选参考星
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@return true=可以继续 RTK 解算；false=本历元跳过（需要时会打印回退提示）
@note 任一系统选不出参考星（RefPrn 为 -1）时返回 false；流动站连单点位置都没有时不打印提示
*/
bool PrepareRelativeObservations(SolverData& rawdata, const SinglePointResult& baseResult,
    const SinglePointResult& roverResult)
{
    // 1、观测数据有效性标记
    MarkObservationsUsable(&rawdata.BasEpk_0, &rawdata.BasEpk_cur);
    MarkObservationsUsable(&rawdata.RovEpk_0, &rawdata.RovEpk_cur);

    // 2、站间单差
    FormInterStationDifference(&rawdata.RovEpk_cur, &rawdata.BasEpk_cur, &rawdata.SdObs);

    // 3、周跳探测
    DetectCycleSlipEvents(&rawdata.SdObs);

    // 4、选参考星，形成双差
    ChooseReferenceSatellite(&rawdata.RovEpk_cur, &rawdata.BasEpk_cur, &rawdata.SdObs,
        &rawdata.DDObs);

    if (rawdata.DDObs.RefPrn[0] == -1 || rawdata.DDObs.RefPrn[1] == -1) {
        // 流动站连单点位置都没有时直接跳过，不打印提示
        if (roverResult.Position[0] == 0 && roverResult.Position[1] == 0 &&
            roverResult.Position[2] == 0)
            return false;

        PrintSppFallback(rawdata.RovEpk_cur.Time, baseResult, roverResult, "");
        return false;
    }
    return true;
}

/*!
@brief 在实时模式下打印一条“退回单点定位”的控制台提示
@param roverTime 当前流动站历元时间
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param pdopNote 跟在 PDOP 数值后面的说明文字，可为空字符串
@return 无
@note 文件模式（ANinfo.IsFileData != 0）下不输出任何内容
*/
void PrintSppFallback(const GpsTime& roverTime, const SinglePointResult& baseResult,
    const SinglePointResult& roverResult, const char* pdopNote)
{
    if (ANinfo.IsFileData != 0) return;                 // 仅实时模式打印

    double dxyz[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++) dxyz[i] = roverResult.Position[i] - baseResult.Position[i];

    printf("[单点解] 第 %u 周 %.1f 秒   dX=%+0.8f   dY=%+0.8f   dZ=%+0.8f   PDOP=%.3f%s   GPS=%d   BDS=%d   σ=%.3f\n",
        roverTime.Week, roverTime.SecOfWeek,
        dxyz[0], dxyz[1], dxyz[2],
        roverResult.PDOP, pdopNote, roverResult.GPSSatNum, roverResult.BDSSatNum,
        roverResult.SigmaPos);
}

/*!
@brief 把单点定位结果换算成东北天偏差后写入 SPP 结果文件
@param rawdata 原始观测数据，提供基站最佳坐标
@param roverResult 流动站单点定位结果
@param sppFile 已打开的单点定位结果文件
@param sppHeader 表头是否已写出的标志
@return 无
@note 首次调用时写出表头；ENU 偏差以基站坐标（BaseBest.Pos）为参考点
*/
void WriteSppResult(SolverData& rawdata, SinglePointResult& roverResult, FILE* sppFile,
    bool& sppHeader)
{
    if (!sppHeader) { fprintf(sppFile, "#Time X Y Z(m)\n"); sppHeader = true; }

    double basePos[3] = { 0.0, 0.0, 0.0 };
    double denu[3] = { 0.0, 0.0, 0.0 };

    GeodeticToEcef(rawdata.BaseBest.Pos, basePos, R_WGS84, F_WGS84);   // 基站的 ECEF 坐标

    GeodeticCoord baseBlh;
    baseBlh.longitude = rawdata.BaseBest.Pos[0];
    baseBlh.latitude = rawdata.BaseBest.Pos[1];
    baseBlh.height = rawdata.BaseBest.Pos[2];

    EcefOffsetToEnu(basePos, roverResult.Position, &baseBlh, denu);

    WriteSinglePointRecord(&rawdata, &roverResult, denu, sppFile);
}
