#pragma once

#include "gnss_core_types.h"
#include "session.h"

/*!
@brief 同步一个历元的基站/流动站观测数据
@param baseFile 基站观测文件（文件模式）
@param roverFile 流动站观测文件（文件模式）
@param baseSocket 基站数据流（实时模式）
@param roverSocket 流动站数据流（实时模式）
@param rawdata 输出的原始观测数据
@return 1=同步成功，0=时间不齐需要跳过，-1=文件结束，-2=数据流失败
*/
int SynchronizeEpoch(FILE* baseFile, FILE* roverFile, SOCKET& baseSocket, SOCKET& roverSocket,
    SolverData& rawdata);

/*!
@brief 分别对基站、流动站做单点定位与测速（卫星位置也在这里算出）
@param rawdata 原始观测数据，函数会更新其中的单点定位状态
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
*/
void SolveStationPositions(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult);

/*!
@brief 完成进入 RTK 解算前的公共步骤：有效性标记、站间单差、周跳探测、选参考星
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@return true=可以继续 RTK 解算；false=本历元跳过（需要时会打印回退提示）
*/
bool PrepareRelativeObservations(SolverData& rawdata, const SinglePointResult& baseResult,
    const SinglePointResult& roverResult);

/*!
@brief 在实时模式下打印一条“退回单点定位”的控制台提示
@param roverTime 当前流动站历元时间
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param pdopNote 跟在 PDOP 数值后面的说明文字，可为空字符串
*/
void PrintSppFallback(const GpsTime& roverTime, const SinglePointResult& baseResult,
    const SinglePointResult& roverResult, const char* pdopNote);

/*!
@brief 把单点定位结果换算成东北天偏差后写入 SPP 结果文件
@param rawdata 原始观测数据，提供基站最佳坐标
@param roverResult 流动站单点定位结果
@param sppFile 已打开的单点定位结果文件
@param sppHeader 表头是否已写出的标志
*/
void WriteSppResult(SolverData& rawdata, SinglePointResult& roverResult, FILE* sppFile,
    bool& sppHeader);

/*!
@brief 处理一个卡尔曼滤波历元
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param ekf 卡尔曼滤波状态（跨历元保留）
@param out 输出文件与表头状态
@param stats 会话统计量
*/
void ProcessEkfEpoch(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult, KalmanState& ekf, SessionOutputs& out, SessionStats& stats);

/*!
@brief 处理一个最小二乘历元
@param rawdata 原始观测数据
@param baseResult 基站单点定位结果
@param roverResult 流动站单点定位结果
@param floatSolution 浮点解结果（跨历元保留）
@param ekf 卡尔曼滤波状态，仅用于绘图格式参数（最小二乘模式不使用其中数值）
@param out 输出文件与表头状态
@param stats 会话统计量
*/
void ProcessLsqEpoch(SolverData& rawdata, SinglePointResult& baseResult,
    SinglePointResult& roverResult, FloatSolution& floatSolution, KalmanState& ekf,
    SessionOutputs& out, SessionStats& stats);
