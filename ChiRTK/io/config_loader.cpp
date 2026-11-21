#include"gnss_core_types.h"


namespace
{
    /*!
    @brief 取冒号之后、去掉首部空白后的配置值文本
    @param line 已截断注释的一行配置文本
    @return 指向配置值起始位置的指针，无冒号时返回 nullptr
    */
    const char* ValueAfterColon(const char* line)
    {
        const char* colon = strchr(line, ':');
        if (colon == NULL) return NULL;

        const char* value = colon + 1;
        while (*value == ' ' || *value == '\t') value++;
        return value;
    }

    /*!
    @brief 报告解析失败
    @param what 配置项的中文名称，用于错误提示
    @param ConfigName 内置配置的名称
    @return 固定返回 false，供调用方直接返回
    */
    bool FailConfigItem(const char* what, const char* ConfigName)
    {
        printf("  配置项“%s”格式错误（内置配置：%s）。\n", what, ConfigName);
        return false;
    }

    /*!
    @brief 解析成功后打印中文配置汇总面板
    @param config 已填好的接收机配置
    */
    void PrintConfigSummary(const RoverConfig& config)
    {
        cout << "  ----------------------------------------------------------------\n";
        cout << "    本次运行配置\n";
        cout << "  ----------------------------------------------------------------\n";
        if (config.IsFileData == 1) {
            cout << "    数据来源：观测文件（事后处理）\n";
            printf("    基站观测文件：%s\n", config.BasObsDatFile);
            printf("    流动站观测文件：%s\n", config.RovObsDatFile);
        }
        else {
            cout << "    数据来源：网络数据流（实时处理）\n";
            printf("    基站地址：%s:%d\n", config.BasNetIP, config.BasNetPort);
            printf("    流动站地址：%s:%d\n", config.RovNetIP, config.RovNetPort);
            printf("    流动站串口：COM%d，波特率 %d\n", config.RovPort, config.RovBaud);
        }
        if (config.RTKProcMode == 1) {
            cout << "    解算方法：卡尔曼滤波（EKF）\n";
            printf("    滤波结果文件：%s\n", config.EkfFile);
        }
        else {
            cout << "    解算方法：最小二乘（LSQ）\n";
            printf("    浮点解结果文件：%s\n", config.FloatFile);
            printf("    固定解结果文件：%s\n", config.FixedFile);
        }
        printf("    单点定位结果文件：%s\n", config.SppFile);
        printf("    观测噪声：伪距 %.3f，载波 %.3f\n", config.CodeNoise, config.CPNoise);
        printf("    高度截止角：%.1f 度\n", config.ElevThreshold);
        printf("    固定解 Ratio 阈值：%.1f\n", config.RatioThres);
        cout << "  ----------------------------------------------------------------\n";
    }
}

/*!
@brief 解析编译进程序的内置配置文本，填充流动站配置结构体的各配置项
@param ConfigName 内置配置的名称（如 LSQ_FILE），仅用于日志和错误提示
@param ConfigText 内置配置文本，格式为"键名: 值 ! 注释"，行间以换行分隔
@param ANInfo 输出的配置信息结构体指针
@return true=解析成功 false=配置项格式非法或文本为空
@note 行内以 "!" 作为注释起始符；文件或通信模式决定哪些配置项被解析；
      解析时先按冒号切出配置值，再交给 sscanf 按类型转换，避免逐条重复写格式串
*/
bool ParseRoverConfig(const char ConfigName[], const char ConfigText[], RoverConfig* ANInfo)
{
    const char CfgStr[CFGITEMS][32] = {
        "FROM FILE OR COM",             // 配置项1：指定是否为文件（FILE）或通信（COM）模式
        "RTK PROCESSING MODE",          // 配置项2：RTK处理模式
        "ROVER IP ADDRESS AND PORT",    // 配置项3：流动站的IP地址和端口
        "BASE IP ADDRESS AND PORT",     // 配置项4：基站的IP地址和端口
        "ROVER COM SETUP",              // 配置项5：流动站的串口设置
        "BASE OBSDATA SOURCE FILE",     // 配置项6：基站的观测数据源文件
        "ROVER OBSDATA SOURCE FILE",    // 配置项7：流动站的观测数据源文件
        "POSITION SPP FILE",            // 配置项8：单点定位结果文件
        "POSITION FLOAT FILE",          // 配置项9：浮点解定位结果文件
        "POSITION FIXED FILE",          // 配置项10：固定解定位结果文件
        "POSITION EKF FILE",            // 配置项11：卡尔曼滤波定位结果文件
        "CODE AND CARRIER PHASE NOISE", // 配置项12：伪距噪声和载波相位噪声
        "THRESHOLD FOR ELEVATION MASK", // 配置项13：卫星高度角的阈值
        "RATIO FOR DD FIXED SOLUTION"   // 配置项14：差分固定解的比值阈值
    };

    if (ConfigText == NULL || ANInfo == NULL)
    {
        printf("  运行失败：内置配置为空（%s）。\n", ConfigName);
        return false;
    }

    char Line[256];

    // 逐行扫描配置文本
    const char* cursor = ConfigText;
    while (*cursor != '\0')
    {
        const char* eol = strchr(cursor, '\n');              // 行尾
        const size_t rawLen = (eol != NULL) ? static_cast<size_t>(eol - cursor) : strlen(cursor);
        size_t copyLen = rawLen;
        if (copyLen >= sizeof(Line)) copyLen = sizeof(Line) - 1;   // 超长行按缓冲区长度截断

        memcpy(Line, cursor, copyLen);
        Line[copyLen] = 0;

        cursor = (eol != NULL) ? eol + 1 : cursor + rawLen;  // 指向下一行

        char* comment = strchr(Line, '!');           // 注释起始位置
        if (comment != NULL) *comment = 0;           // 截掉注释部分

        const char* value = NULL;                    // 该行的配置值
        int item = -1;                               // 命中的配置项下标

        for (int i = 0; i < CFGITEMS; i++)
        {
            if (strstr(Line, CfgStr[i]) != NULL)
            {
                item = i;
                value = ValueAfterColon(Line);
                break;                               // 一行只对应一个配置项
            }
        }

        if (item < 0 || value == NULL) continue;     // 空行、纯注释行或未知行

        switch (item)
        {
        case 0:                                      // 文件 / 通信模式
            if (sscanf(value, "%d", &ANInfo->IsFileData) != 1)
                return FailConfigItem("数据来源（0=网络，1=文件）", ConfigName);
            break;

        case 1:                                      // RTK 处理模式
            if (sscanf(value, "%d", &ANInfo->RTKProcMode) != 1)
                return FailConfigItem("解算模式（1=EKF，2=LSQ）", ConfigName);
            break;

        case 2:                                      // 流动站 IP 与端口（仅通信模式）
            if (ANInfo->IsFileData != 0) break;
            if (sscanf(value, "%s %d", ANInfo->RovNetIP, &ANInfo->RovNetPort) != 2)
                return FailConfigItem("流动站 IP 与端口", ConfigName);
            break;

        case 3:                                      // 基站 IP 与端口（仅通信模式）
            if (ANInfo->IsFileData != 0) break;
            if (sscanf(value, "%s %d", ANInfo->BasNetIP, &ANInfo->BasNetPort) != 2)
                return FailConfigItem("基站 IP 与端口", ConfigName);
            break;

        case 4:                                      // 流动站串口设置（仅通信模式）
            if (ANInfo->IsFileData != 0) break;
            if (sscanf(value, "%d %d", &ANInfo->RovPort, &ANInfo->RovBaud) != 2)
                return FailConfigItem("流动站串口", ConfigName);
            break;

        case 5:                                      // 基站观测数据文件（仅文件模式）
            if (ANInfo->IsFileData != 1) break;
            if (sscanf(value, "%s", ANInfo->BasObsDatFile) != 1)
                return FailConfigItem("基站观测文件", ConfigName);
            break;

        case 6:                                      // 流动站观测数据文件（仅文件模式）
            if (ANInfo->IsFileData != 1) break;
            if (sscanf(value, "%s", ANInfo->RovObsDatFile) != 1)
                return FailConfigItem("流动站观测文件", ConfigName);
            break;

        case 7:                                      // SPP 结果文件
            if (sscanf(value, "%s", ANInfo->SppFile) != 1)
                return FailConfigItem("单点定位结果文件", ConfigName);
            break;

        case 8:                                      // 浮点解结果文件（仅 LSQ 模式）
            if (ANInfo->RTKProcMode != 2) break;
            if (sscanf(value, "%s", ANInfo->FloatFile) != 1)
                return FailConfigItem("浮点解结果文件", ConfigName);
            break;

        case 9:                                      // 固定解结果文件（仅 LSQ 模式）
            if (ANInfo->RTKProcMode != 2) break;
            if (sscanf(value, "%s", ANInfo->FixedFile) != 1)
                return FailConfigItem("固定解结果文件", ConfigName);
            break;

        case 10:                                     // EKF 结果文件（仅 EKF 模式）
            if (ANInfo->RTKProcMode != 1) break;
            if (sscanf(value, "%s", ANInfo->EkfFile) != 1)
                return FailConfigItem("滤波结果文件", ConfigName);
            break;

        case 11:                                     // 观测值噪声
            if (sscanf(value, "%lf %lf", &ANInfo->CodeNoise, &ANInfo->CPNoise) != 2)
                return FailConfigItem("伪距与载波噪声", ConfigName);
            break;

        case 12:                                     // 高度角阈值
            if (sscanf(value, "%lf", &ANInfo->ElevThreshold) != 1)
                return FailConfigItem("高度截止角", ConfigName);
            break;

        case 13:                                     // Ratio 检验阈值
            if (sscanf(value, "%lf", &ANInfo->RatioThres) != 1)
                return FailConfigItem("固定解 Ratio 阈值", ConfigName);
            break;

        default:
            break;
        }
    }

    cout << "\n";
    PrintConfigSummary(*ANInfo);

    return true;
}
