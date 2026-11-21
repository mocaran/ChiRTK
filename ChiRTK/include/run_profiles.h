#pragma once

#include "gnss_core_types.h"

// The menu and its built-in configurations live in one place.
struct RunProfile {
    const char* configName;     // 内置配置的名称，仅用于日志与错误提示
    const char* configText;     // 编译进程序的内置配置文本，见 builtin_config.h
    const char* description;
    const char* resultDir;
};

/*!
@brief 按下标从内置运行模式表中取出对应的内置配置与结果目录
@param choice 运行模式编号，有效范围0~3（0=LSQ文件，1=LSQ实时，2=EKF文件，3=EKF实时）
@return 指向对应RunProfile常量的指针；choice越界时返回nullptr
@note 返回的指针指向静态常量表，调用方不应释放或修改
*/
const RunProfile* LookupRunProfile(int choice);

// Parse the selected built-in configuration and check inputs before creating output files.
/*!
@brief 解析指定模式的内置配置并检查输入文件，创建结果目录
@param mode 运行模式，提供内置配置文本与结果目录
@param config 输出的接收机配置信息，解析成功后填满
@return true=准备成功 false=配置解析失败或输入文件缺失
@note 文件数据模式下要求基站与流动站观测文件均已存在
*/
bool PrepareRunProfile(const RunProfile& mode, RoverConfig& config);
