#include "run_profiles.h"

#include "builtin_config.h"

#include <array>
#include <filesystem>
#include <iostream>
#include <system_error>

namespace {
// 四个运行模式与四份内置配置一一对应，配置文本见 builtin_config.h
// 编号 1、2 使用观测文件，编号 3、4 使用网络数据流；1/3 为最小二乘，2/4 为卡尔曼滤波
constexpr std::array<RunProfile, 4> kModes{ {
    {"LSQ_FILE", kBuiltinConfigLsqFile, "观测文件 + 最小二乘（LSQ），事后处理", "Result\\LSQ_FILE"},
    {"EKF_FILE", kBuiltinConfigEkfFile, "观测文件 + 卡尔曼滤波（EKF），事后处理", "Result\\EKF_FILE"},
    {"LSQ_LIVE", kBuiltinConfigLsqLive, "网络数据流 + 最小二乘（LSQ），实时处理", "Result\\LSQ_LIVE"},
    {"EKF_LIVE", kBuiltinConfigEkfLive, "网络数据流 + 卡尔曼滤波（EKF），实时处理", "Result\\EKF_LIVE"},
} };

// 文件数据模式下必须存在的输入，按"名称 + 路径"成对记录，便于统一检查
struct RequiredInput {
    const char* label;
    const char* path;
};
}

/*!
@brief 按下标从内置运行模式表中取出对应的内置配置与结果目录
@param choice 运行模式编号，有效范围1~4（1=LSQ文件，2=EKF文件，3=LSQ实时，4=EKF实时）
@return 指向对应RunProfile常量的指针；choice越界时返回nullptr
@note 返回的指针指向静态常量表，调用方不应释放或修改
*/
const RunProfile* LookupRunProfile(int choice)
{
    if (choice < 1 || choice > static_cast<int>(kModes.size())) return nullptr;
    return &kModes[static_cast<size_t>(choice - 1)];
}

/*!
@brief 解析指定模式的内置配置并检查输入文件，创建结果目录
@param mode 运行模式，提供内置配置文本与结果目录
@param config 输出的接收机配置信息，解析成功后填满
@return true=准备成功，false=配置解析失败或输入文件缺失
@note 文件数据模式下要求基站与流动站观测文件均已存在；
      结果目录用带 error_code 的重载创建，失败时直接取回错误描述而不用异常
*/
bool PrepareRunProfile(const RunProfile& mode, RoverConfig& config)
{
    if (!ParseRoverConfig(mode.configName, mode.configText, &config)) {
        std::cerr << "  运行失败：内置配置解析未通过：" << mode.configName << '\n';
        return false;
    }

    if (config.IsFileData == 1) {
        const RequiredInput inputs[] = {
            {"基站", config.BasObsDatFile},
            {"流动站", config.RovObsDatFile},
        };

        bool allPresent = true;

        for (const RequiredInput& input : inputs) {
            if (!std::filesystem::exists(input.path)) {
                if (allPresent) std::cerr << "  运行失败：缺少观测文件：\n";
                std::cerr << "    " << input.label << "：" << input.path << '\n';
                allPresent = false;
            }
        }

        if (!allPresent) return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(mode.resultDir, ec);

    if (ec) {
        std::cerr << "  运行失败：无法创建结果目录：" << ec.message() << '\n';
        return false;
    }
    return true;
}
