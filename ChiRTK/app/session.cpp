#include "session.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

#include "epoch_processing.h"



namespace {

using Clock = std::chrono::steady_clock;

// 文件模式下保存输入文件句柄；实时模式下保存两条数据流
struct SessionFiles
{
    FILE* baseFile = nullptr;
    FILE* roverFile = nullptr;
    SOCKET baseSocket = INVALID_SOCKET;
    SOCKET roverSocket = INVALID_SOCKET;
};

struct CloseStatus
{
    bool readFailed = false;
    bool writeFailed = false;
};

enum class LoopResult
{
    SourceEnded,
    TimeReached,
    StreamFailed
};

/*!
@brief 按“追加/覆盖”策略打开一个结果文件，并判断它是否还没有表头
@param path 文件路径
@param append true=追加写（实时模式），false=覆盖写（文件模式）
@param isNew 输出：true 表示文件是新建或空的，需要写表头
@return 文件句柄，失败时为 nullptr
*/
FILE* OpenResultFile(const std::string& path, bool append, bool& isNew)
{
    std::error_code ec;
    isNew = true;
    if (append) {
        isNew = !std::filesystem::exists(path, ec) || std::filesystem::file_size(path, ec) == 0;
    }
    return fopen(path.c_str(), append ? "a" : "w");
}

/*!
@brief 打开本次会话需要的全部输入/输出文件，并逐一检查
@param mode 运行模式，提供结果目录
@param out 输出文件与表头状态
@param files 输入文件句柄
@return true=全部就绪，false=有文件打开失败（已打印原因）
*/
bool OpenSessionFiles(const RunProfile& mode, SessionOutputs& out, SessionFiles& files)
{
    const bool appendResults = (ANinfo.IsFileData == 0);
    const std::string plotPath = std::string(mode.resultDir) + "\\RTKPlot.pos";
    std::cout << "\n    结果目录：" << mode.resultDir
        << (appendResults ? "（新结果追加写入）\n" : "（同名旧文件将被覆盖）\n");

    bool sppNew = true, floatNew = true, fixedNew = true, ekfNew = true, plotNew = true;

    // 打开顺序与原实现保持一致：绘图文件、观测文件、SPP、浮点/固定、卡尔曼
    if (ANinfo.RTKProcMode == 1) out.ekfPlot = OpenResultFile(plotPath, appendResults, plotNew);
    if (ANinfo.RTKProcMode == 2) out.lsqPlot = OpenResultFile(plotPath, appendResults, plotNew);
    files.baseFile = ANinfo.IsFileData == 1 ? fopen(ANinfo.BasObsDatFile, "rb") : nullptr;
    files.roverFile = ANinfo.IsFileData == 1 ? fopen(ANinfo.RovObsDatFile, "rb") : nullptr;
    out.spp = OpenResultFile(ANinfo.SppFile, appendResults, sppNew);
    if (ANinfo.RTKProcMode == 2) out.floatSolution = OpenResultFile(ANinfo.FloatFile, appendResults, floatNew);
    if (ANinfo.RTKProcMode == 2) out.fixedSolution = OpenResultFile(ANinfo.FixedFile, appendResults, fixedNew);
    if (ANinfo.RTKProcMode == 1) out.ekf = OpenResultFile(ANinfo.EkfFile, appendResults, ekfNew);

    if (ANinfo.IsFileData == 1 && files.baseFile == nullptr) {
        fprintf(stderr, "  运行失败：无法打开基站观测文件。\n");
        return false;
    }
    if (ANinfo.IsFileData == 1 && files.roverFile == nullptr) {
        fprintf(stderr, "  运行失败：无法打开流动站观测文件。\n");
        return false;
    }
    if (out.spp == nullptr) {
        fprintf(stderr, "  运行失败：无法创建单点定位结果文件。\n");
        return false;
    }
    if (ANinfo.RTKProcMode == 2) {
        if (out.floatSolution == nullptr) {
            fprintf(stderr, "  运行失败：无法创建浮点解结果文件。\n");
            return false;
        }
        if (out.fixedSolution == nullptr) {
            fprintf(stderr, "  运行失败：无法创建固定解结果文件。\n");
            return false;
        }
        if (out.lsqPlot == nullptr) {
            fprintf(stderr, "  运行失败：无法创建最小二乘绘图文件。\n");
            return false;
        }
    }
    if (ANinfo.RTKProcMode == 1 && out.ekf == nullptr) {
        fprintf(stderr, "  运行失败：无法创建卡尔曼滤波结果文件。\n");
        return false;
    }
    if (ANinfo.RTKProcMode == 1 && out.ekfPlot == nullptr) {
        fprintf(stderr, "  运行失败：无法创建卡尔曼滤波绘图文件。\n");
        return false;
    }

    // 追加模式下 true 表示文件里已有表头，不再重复写
    out.sppHeader = !sppNew;
    out.floatHeader = !floatNew;
    out.fixedHeader = !fixedNew;
    out.ekfHeader = !ekfNew;
    return true;
}

/*!
@brief 实时模式：询问采集时长并连接基站、流动站两条数据流
@param files 数据流句柄
@param realtimeMinutes 输出：用户输入的采集分钟数
@return true=连接成功，false=输入或连接失败（已打印原因）
*/
bool PrepareLiveStreams(SessionFiles& files, double& realtimeMinutes)
{
    std::cout << "\n    请输入实时采集时长（分钟）> ";
    if (!(std::cin >> realtimeMinutes) || realtimeMinutes <= 0.0) {
        std::cerr << "  运行失败：采集时长必须大于 0。\n";
        return false;
    }

    std::cout << "    正在连接基站数据流 " << ANinfo.BasNetIP << ':' << ANinfo.BasNetPort << " ...\n";
    if (!ConnectStreamSocket(files.baseSocket, ANinfo.BasNetIP, ANinfo.BasNetPort)) {
        std::cerr << "  运行失败：无法连接基站数据流。\n";
        return false;
    }

    std::cout << "    正在连接流动站数据流 " << ANinfo.RovNetIP << ':' << ANinfo.RovNetPort << " ...\n";
    if (!ConnectStreamSocket(files.roverSocket, ANinfo.RovNetIP, ANinfo.RovNetPort)) {
        DisconnectStreamSocket(files.baseSocket);
        std::cerr << "  运行失败：无法连接流动站数据流。\n";
        return false;
    }
    return true;
}

/*!
@brief 每 5 秒打印一次运行进度；文件模式显示百分比，实时模式显示已用时间
*/
class ProgressReporter
{
public:
    /*!
    @brief 保存进度打印所需的输入文件与采集时长信息，并记录开始时刻
    @param fileMode true=文件模式
    @param roverFile 流动站观测文件（文件模式），用于按字节数估算进度
    @param roverFileBytes 流动站观测文件总字节数
    @param sizeKnown 文件大小是否已知
    @param realtimeMinutes 实时模式采集分钟数
    */
    ProgressReporter(bool fileMode, FILE* roverFile, std::uintmax_t roverFileBytes,
        bool sizeKnown, double realtimeMinutes)
        : fileMode_(fileMode), roverFile_(roverFile), roverFileBytes_(roverFileBytes),
          sizeKnown_(sizeKnown), realtimeMinutes_(realtimeMinutes),
          startedAt_(Clock::now()), lastProgress_(startedAt_)
    {
    }

    /*!
    @brief 按 5 秒节拍打印一行运行进度；文件模式显示百分比，实时模式显示已用时间
    @param stats 当前会话统计量
    */
    void Update(const SessionStats& stats)
    {
        const auto now = Clock::now();
        if (stats.syncedEpochs != 1 && now - lastProgress_ < std::chrono::seconds(5)) return;

        const auto seconds =
            std::chrono::duration_cast<std::chrono::seconds>(now - startedAt_).count();
        if (fileMode_ && sizeKnown_ && roverFileBytes_ > 0) {
            const long bytesRead = ftell(roverFile_);
            double percent = 0.0;
            if (bytesRead >= 0) {
                const std::uintmax_t readBytes = static_cast<std::uintmax_t>(bytesRead);
                percent = readBytes >= roverFileBytes_ ? 100.0
                    : 100.0 * static_cast<double>(readBytes) / static_cast<double>(roverFileBytes_);
            }
            printf("    进度 %5.1f%%   已用 %lld 秒   同步 %llu 个历元   已保存 %llu 个结果\n",
                percent, static_cast<long long>(seconds), stats.syncedEpochs, stats.savedEpochs);
        }
        else {
            printf("    进度   已用 %lld 秒 / 共 %.0f 秒   同步 %llu 个历元   已保存 %llu 个结果\n",
                static_cast<long long>(seconds), realtimeMinutes_ * 60.0,
                stats.syncedEpochs, stats.savedEpochs);
        }
        fflush(stdout);
        lastProgress_ = now;
    }

    /*!
    @brief 返回自会话开始以来经过的秒数
    */
    double ElapsedSeconds() const
    {
        return std::chrono::duration<double>(Clock::now() - startedAt_).count();
    }

private:
    bool fileMode_;
    FILE* roverFile_;
    std::uintmax_t roverFileBytes_;
    bool sizeKnown_;
    double realtimeMinutes_;
    Clock::time_point startedAt_;
    Clock::time_point lastProgress_;
};

/*!
@brief 逐历元主循环：取数、单点定位、RTK 解算，直到文件结束或采集时间到
@param out 输出文件与表头状态
@param files 输入文件/数据流
@param stats 会话统计量
@param progress 进度打印器
@param fileMode true=文件模式
@param realtimeMinutes 实时模式采集分钟数
@return 循环结束原因
*/
LoopResult RunEpochLoop(SessionOutputs& out, SessionFiles& files, SessionStats& stats,
    ProgressReporter& progress, bool fileMode, double realtimeMinutes)
{
    // 这些状态要跨历元保留：上一历元观测、滤波状态、模糊度信息等
    SolverData rawdata;
    SinglePointResult baseResult;
    SinglePointResult roverResult;
    FloatSolution floatSolution;
    KalmanState ekf;

    while (true) {
        if (!fileMode && progress.ElapsedSeconds() >= realtimeMinutes * 60.0)
            return LoopResult::TimeReached;

        // 1、时间同步
        const int syncFlag = SynchronizeEpoch(files.baseFile, files.roverFile,
            files.baseSocket, files.roverSocket, rawdata);
        if (syncFlag == 0) { progress.Update(stats); continue; }   // 两侧时间不齐，跳过
        if (syncFlag == -1) return LoopResult::SourceEnded;
        if (syncFlag == -2) return LoopResult::StreamFailed;
        ++stats.syncedEpochs;
        progress.Update(stats);

        // 2、基站、流动站单点定位
        SolveStationPositions(rawdata, baseResult, roverResult);

        // 3、单双差准备；不能进入 RTK 的历元直接跳过
        if (!PrepareRelativeObservations(rawdata, baseResult, roverResult)) continue;

        // 4、按配置走卡尔曼滤波或最小二乘
        if (ANinfo.RTKProcMode == 1) {
            ProcessEkfEpoch(rawdata, baseResult, roverResult, ekf, out, stats);
        }
        else {
            ProcessLsqEpoch(rawdata, baseResult, roverResult, floatSolution, ekf, out, stats);
        }
    }
}

/*!
@brief 收尾：检查输入读错误，关闭观测文件/数据流和结果文件
@param out 输出文件
@param files 输入文件/数据流
@param fileMode true=文件模式
@return 读、写是否出错
*/
CloseStatus CloseSession(SessionOutputs& out, SessionFiles& files, bool fileMode)
{
    CloseStatus status;
    if (fileMode) {
        status.readFailed = ferror(files.baseFile) != 0 || ferror(files.roverFile) != 0;
        fclose(files.baseFile);
        fclose(files.roverFile);
    }
    else {
        DisconnectStreamSocket(files.baseSocket);
        DisconnectStreamSocket(files.roverSocket);
    }

    status.writeFailed |= fclose(out.spp) != 0;
    if (ANinfo.RTKProcMode == 2) {
        status.writeFailed |= fclose(out.floatSolution) != 0;
        status.writeFailed |= fclose(out.fixedSolution) != 0;
        status.writeFailed |= fclose(out.lsqPlot) != 0;
    }
    else {
        status.writeFailed |= fclose(out.ekf) != 0;
        status.writeFailed |= fclose(out.ekfPlot) != 0;
    }
    return status;
}

/*!
@brief 打印最终结果，并换算成程序退出码
@param mode 运行模式，提供结果目录
@param stats 会话统计量
@param result 循环结束原因
@param close 读、写错误状态
@param elapsed 总耗时（秒）
@return 0=成功，1=失败
*/
int ReportSession(const RunProfile& mode, const SessionStats& stats, LoopResult result,
    const CloseStatus& close, double elapsed)
{
    if (result == LoopResult::StreamFailed || close.readFailed || close.writeFailed ||
        stats.fixedEpochs + stats.floatEpochs == 0) {
        const char* reason = nullptr;
        if (result == LoopResult::StreamFailed) {
            reason = "网络数据流中断或返回错误";
        }
        else if (close.readFailed) {
            reason = "观测文件未能完整读取";
        }
        else if (close.writeFailed) {
            reason = "结果文件未能完整写入";
        }
        else if (stats.savedEpochs == 0) {
            reason = "没有保存任何定位结果，请检查基站与流动站的观测是否重叠";
        }
        else {
            reason = "只保存了单点定位结果，未得到 RTK 浮点解或固定解";
        }

        printf("\n  ================================================================\n");
        printf("    解算失败\n");
        printf("  ----------------------------------------------------------------\n");
        printf("    失败原因：%s\n", reason);
        printf("    时间同步历元：%llu\n", stats.syncedEpochs);
        printf("    保存结果个数：%llu\n", stats.savedEpochs);
        printf("    总耗时：%.1f 秒\n", elapsed);
        printf("    结果目录中可能存在不完整或空文件：%s\n", mode.resultDir);
        printf("  ================================================================\n");
        return 1;
    }

    const char* endReason = "处理被提前终止";
    if (result == LoopResult::SourceEnded) {
        endReason = "观测文件已读取完毕";
    }
    else if (result == LoopResult::TimeReached) {
        endReason = "已达到设定的实时采集时长";
    }

    printf("\n  ================================================================\n");
    printf("    解算成功\n");
    printf("  ----------------------------------------------------------------\n");
    printf("    结束原因：%s\n", endReason);
    printf("    时间同步历元：%llu\n", stats.syncedEpochs);
    printf("    保存结果个数：%llu（固定解 %llu / 浮点解 %llu / 单点解 %llu）\n",
        stats.savedEpochs, stats.fixedEpochs, stats.floatEpochs, stats.sppEpochs);
    printf("    跳过历元个数：%llu\n", stats.syncedEpochs - stats.savedEpochs);
    printf("    总耗时：%.1f 秒\n", elapsed);
    printf("    结果目录：%s\n", mode.resultDir);
    printf("  ================================================================\n");

    // 机器可读统计行：analysis\scripts\analyze_run_stats.py 按此格式解析，键名保持英文
    printf("RUN SUCCEEDED: synchronized=%llu, saved=%llu (fixed=%llu, float=%llu, SPP=%llu), skipped=%llu, elapsed=%.1fs.\n",
        stats.syncedEpochs, stats.savedEpochs, stats.fixedEpochs, stats.floatEpochs,
        stats.sppEpochs, stats.syncedEpochs - stats.savedEpochs, elapsed);
    return 0;
}

}  // namespace

/*!
@brief 运行一次完整的解算会话（打开文件、逐历元处理、收尾报告）
@param mode 已选中的运行模式，提供内置配置与结果目录
@return 程序退出码：0=成功，1=输入、数据流或输出失败
@note 内部按“打开输出文件 → 建立实时连接 → 逐历元循环 → 关闭文件 → 汇总报告”的顺序执行
*/
int RunSession(const RunProfile& mode)
{
    SessionOutputs outputs;
    SessionFiles files;
    if (!OpenSessionFiles(mode, outputs, files)) return 1;

    double realtimeMinutes = 0.0;
    if (ANinfo.IsFileData == 0 && !PrepareLiveStreams(files, realtimeMinutes)) return 1;

    const bool fileMode = (ANinfo.IsFileData == 1);
    std::error_code sizeError;
    const std::uintmax_t roverFileBytes = fileMode
        ? std::filesystem::file_size(ANinfo.RovObsDatFile, sizeError) : 0;

    ProgressReporter progress(fileMode, files.roverFile, roverFileBytes, !sizeError,
        realtimeMinutes);
    std::cout << "\n    开始逐历元解算，运行期间无需任何输入。\n" << std::flush;

    SessionStats stats;
    const LoopResult result = RunEpochLoop(outputs, files, stats, progress, fileMode,
        realtimeMinutes);
    const CloseStatus close = CloseSession(outputs, files, fileMode);
    return ReportSession(mode, stats, result, close, progress.ElapsedSeconds());
}
