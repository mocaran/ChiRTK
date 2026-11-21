#include "gnss_core_types.h"
#include "run_profiles.h"
#include "session.h"

#include <iostream>

RoverConfig ANinfo;   // 当前运行配置，解算模块通过 session.h 中的 extern 声明访问
int codemode = -1;    // 1=文件模式，0=实时模式

namespace {

/*!
@brief 打印运行模式选择面板（1、2 为观测文件，3、4 为网络数据流）
*/
void PrintModeMenu()
{
    std::cout
        << "\n"
        << "  ================================================================\n"
        << "                  ChiRTK  高精度相对定位解算\n"
        << "  ================================================================\n"
        << "    运行模式（输入编号 1-4，按回车确认）\n"
        << "\n"
        << "      1   观测文件   + 最小二乘 LSQ      事后处理\n"
        << "      2   观测文件   + 卡尔曼滤波 EKF    事后处理\n"
        << "      3   网络数据流 + 最小二乘 LSQ      实时处理\n"
        << "      4   网络数据流 + 卡尔曼滤波 EKF    实时处理\n"
        << "\n"
        << "  ----------------------------------------------------------------\n"
        << "    请输入编号 > ";
}

/*!
@brief 从控制台读取运行模式编号
@param choice 输出：用户输入的编号
@return true=读取成功，false=输入结束或不是数字（已打印原因）
*/
bool ReadRunMode(int& choice)
{
    PrintModeMenu();
    if (!(std::cin >> choice)) {
        std::cerr << "\n  运行失败：编号必须是数字（1-4）。\n";
        return false;
    }
    return true;
}

}  // namespace

/*!
@brief 程序入口：选择模式、解析内置配置，然后把逐历元解算交给 RunSession
@return 0=成功，1=输入、配置、数据流或输出失败
*/
int main()
{
    int choice = -1;
    if (!ReadRunMode(choice)) return 1;

    const RunProfile* mode = LookupRunProfile(choice);
    if (mode == nullptr) {
        std::cerr << "\n  运行失败：编号只能取 1、2、3、4，当前输入为 " << choice << "。\n";
        return 1;
    }
    std::cout << "\n    已选择：" << mode->description << "\n";

    if (!PrepareRunProfile(*mode, ANinfo)) {
        std::cerr << "\n  运行失败：内置配置或输入文件检查未通过。\n";
        return 1;
    }
    codemode = ANinfo.IsFileData;

    return RunSession(*mode);
}
