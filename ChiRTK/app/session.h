#pragma once

#include <cstdio>

#include "gnss_core_types.h"
#include "run_profiles.h"

// 主程序里定义的运行配置；各解算模块通过它判断文件/实时模式和结果文件名。
extern RoverConfig ANinfo;

/*!
@brief 一次解算会话中所有输出文件的句柄与表头状态
@note 表头标志为 true 表示文件里已经有表头（实时模式追加写时不再重复写）
*/
struct SessionOutputs
{
    FILE* spp = nullptr;             // 单点定位结果
    FILE* floatSolution = nullptr;   // 浮点解结果
    FILE* fixedSolution = nullptr;   // 固定解结果
    FILE* ekf = nullptr;             // 卡尔曼滤波结果
    FILE* lsqPlot = nullptr;         // 最小二乘绘图数据
    FILE* ekfPlot = nullptr;         // 卡尔曼绘图数据
    bool sppHeader = false;
    bool floatHeader = false;
    bool fixedHeader = false;
    bool ekfHeader = false;
};

/*! @brief 整个会话的逐历元统计量 */
struct SessionStats
{
    unsigned long long syncedEpochs = 0;   // 时间同步成功的历元数
    unsigned long long savedEpochs = 0;    // 实际写出结果的历元数
    unsigned long long fixedEpochs = 0;    // 固定解历元数
    unsigned long long floatEpochs = 0;    // 浮点解历元数
    unsigned long long sppEpochs = 0;      // 只能给出单点定位的历元数
};

/*!
@brief 运行一次完整的解算会话（打开文件、逐历元处理、收尾报告）
@param mode 已选中的运行模式，提供内置配置与结果目录
@return 程序退出码：0=成功，1=输入、数据流或输出失败
*/
int RunSession(const RunProfile& mode);
