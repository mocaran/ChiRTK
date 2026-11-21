#include"gnss_core_types.h"


namespace
{
    /*!
    @brief 按目标类型的字节宽度从内存中取出一个标量
    @tparam T 目标标量类型（double/float/int/unsigned int/short/unsigned short）
    @param p 指向待读取内存的字节指针
    @return 读取到的标量值
    @note 统一用 memcpy 搬运，避免非对齐访问与类型别名问题
    */
    template <typename T>
    inline T ReadScalar(const unsigned char* p)
    {
        T value;
        memcpy(&value, p, sizeof(T));
        return value;
    }

}

/*!
@brief 从指定内存地址读取一个双精度浮点数
@param p 指向待读取内存的字节指针
@return 读取到的双精度浮点数值
@note 使用 memcpy 按 8 字节拷贝，避免非对齐访问与类型别名问题
*/
double ReadAsDouble(unsigned char* p) // 提取 8 字节 double
{
    return ReadScalar<double>(p);
}
/*!
@brief 从指定内存地址读取一个单精度浮点数
@param p 指向待读取内存的字节指针
@return 读取到的单精度浮点数值
@note 使用 memcpy 按 4 字节拷贝，避免非对齐访问与类型别名问题
*/
float ReadAsFloat(unsigned char* p)  // 提取 4 字节 float
{
    return ReadScalar<float>(p);
}
/*!
@brief 从指定内存地址读取一个有符号 32 位整数
@param p 指向待读取内存的字节指针
@return 读取到的有符号 32 位整数值
@note 使用 memcpy 按 4 字节拷贝，避免非对齐访问
*/
int ReadAsInt32(unsigned char* p)    // 提取 4 字节有符号整数
{
    return ReadScalar<int>(p);
}
/*!
@brief 从指定内存地址读取一个无符号 32 位整数
@param p 指向待读取内存的字节指针
@return 读取到的无符号 32 位整数值
@note 使用 memcpy 按 4 字节拷贝，避免非对齐访问
*/
unsigned int ReadAsUint32(unsigned char* p)  // 提取 4 字节无符号整数
{
    return ReadScalar<unsigned int>(p);
}
/*!
@brief 从指定内存地址读取一个有符号 16 位整数
@param p 指向待读取内存的字节指针
@return 读取到的有符号 16 位整数值
@note 使用 memcpy 按 2 字节拷贝，避免非对齐访问
*/
short ReadAsInt16(unsigned char* p)  // 提取 2 字节有符号整数
{
    return ReadScalar<short>(p);
}
/*!
@brief 从指定内存地址读取一个无符号 16 位整数
@param p 指向待读取内存的字节指针
@return 读取到的无符号 16 位整数值
@note 使用 memcpy 按 2 字节拷贝，避免非对齐访问
*/
unsigned short ReadAsUint16(unsigned char* p)    // 提取 2 字节无符号整数
{
    return ReadScalar<unsigned short>(p);
}


/*!
@brief 按位逐比特计算 NovAtel OEM7 二进制消息的 CRC32 校验值
@param buff 待校验数据的首地址指针
@param len 参与校验的字节数，含 28 字节消息头与消息体
@return 计算得到的 32 位 CRC32 校验值
@note 多项式为 POLYCRC32，按位右移逐比特计算，用于消息完整性校验
*/
unsigned int ComputeCrc32(const unsigned char* buff, int len)
{
    int i, j;
    unsigned int crc = 0;
    for (i = 0; i < len; i++)
    {
        crc ^= buff[i];
        for (j = 0; j < 8; j++)
        {
            crc = (crc & 1) ? ((crc >> 1) ^ POLYCRC32) : (crc >> 1);
        }
    }
    return crc;
}


/*!
@brief 把卫星系统与信号类型映射为频率序号和载波波长
@param satSystem 卫星系统标志，0=GPS，4=BDS
@param sigType 信号类型标志
@param freq 输出的频率序号，0=第一频率，1=第二频率
@param wl 输出的载波波长，单位米
@return true=属于所需双频信号，false=应跳过该观测值
@note GPS 取 L1 C/A 与 L2P(Y)，BDS 取 B1I 与 B3I
*/
static bool MapSignalToFrequency(int satSystem, int sigType, int& freq, double& wl)
{
    if (satSystem == 0)                       // GPS
    {
        if (sigType == 0) { freq = 0; wl = WL1_GPS; return true; }
        if (sigType == 9) { freq = 1; wl = WL2_GPS; return true; }
        return false;
    }
    if (satSystem == 4)                       // BDS
    {
        if (sigType == 0 || sigType == 4) { freq = 0; wl = WL1_BDS; return true; }
        if (sigType == 2 || sigType == 6) { freq = 1; wl = WL3_BDS; return true; }
        return false;
    }
    return false;
}

/*!
@brief 解析 RANGE 消息中的伪距、载波相位、多普勒与载噪比观测值
@param data 指向完整消息起始位置的字节缓冲区
@param obs 输出的当前历元观测值结构体指针
@return 无
@note 消息体每条观测值固定 44 字节；仅保留 GPS 与 BDS 的 L1/L2、B1I/B3I 双频信号
*/
void ParseRangeMessage(unsigned char* data, EpochObservation* obs)
{
    unsigned char* p = data + 28;   // 跳过同步头和消息头，直接指向消息体

    //1、从消息头中解码得到观测时刻，该时刻为接收机钟表面时，用GpsTime结构体表示。
    obs->Time.Week = ReadAsUint16(data + 14);
    obs->Time.SecOfWeek = ReadAsUint32(data + 16) * 1E-3;    // 毫秒要乘e-3

    //2、解码得到观测值数量，为所有卫星所有信号观测值的总数。
    const int obsNum = (int)ReadAsUint32(p);
    memset(obs->SatObs, 0, MAXCHANNUM * sizeof(SatelliteObservation));    // 将当前历元的所有卫星观测数据清空

    int lastAllocated = 0;    // 最近一次分配的空槽下标
    p += 4;                   // 跳过 ObsNum 这 4 个字节

    //3、对所有信号观测值进行循环解码
    for (int i = 0; i < obsNum; i++, p += 44)
    {
        //① 从跟踪状态标记中取出各标志位
        const unsigned int chanStatus = ReadAsUint32(p + 40);
        const int parityFlag = (chanStatus >> 11) & 0x01;
        const int phaseLockFlag = (chanStatus >> 10) & 0x01;
        const int codeLockedFlag = (chanStatus >> 12) & 0x01;
        const int satSystem = (chanStatus >> 16) & 0x07;
        const int sigType = (chanStatus >> 21) & 0x1F;

        //② ③ 系统/信号类型筛选，并记录频率序号与波长
        int freq = 0;
        double wl = 0.0;
        if (!MapSignalToFrequency(satSystem, sigType, freq, wl)) continue;

        GnssSystem sys = (satSystem == 0) ? GPS : BDS;

        //④ 在当前观测值结构体中查找同一卫星的槽位，找不到则占用一个空槽
        const int prn = ReadAsUint16(p);
        int slot = -1;

        for (int j = 0; j < MAXCHANNUM; j++)
        {
            if (obs->SatObs[j].System == sys && obs->SatObs[j].Prn == prn)
            {
                slot = j;
                break;
            }
            if (obs->SatObs[j].Prn == 0)
            {
                slot = j;
                lastAllocated = j;
                break;
            }
        }

        if (slot < 0) continue;    // 槽位已满，丢弃该观测值

        SatelliteObservation* sat = &obs->SatObs[slot];

        sat->Prn = prn;                       // 卫星编号
        sat->System = sys;                    // 所属系统

        sat->P[freq] = (codeLockedFlag == 1) ? ReadAsDouble(p + 4) : 0.0;          // 伪距P
        sat->L[freq] = -wl * ((phaseLockFlag == 1) ? ReadAsDouble(p + 16) : 0.0);  // 载波相位L
        sat->D[freq] = -wl * ReadAsFloat(p + 28);                                  // 多普勒频率D
        sat->cn0[freq] = ReadAsFloat(p + 32);                                      // 载噪比
        sat->LockTime[freq] = ReadAsFloat(p + 36);                                 // 锁定时间
        sat->half[freq] = parityFlag;                                                                          // 半周模糊度标志
    }

    //4. 统计有效的卫星数量
    obs->SatNum = lastAllocated + 1;
}


/*!
@brief 解析 NovAtel OEM7 GPS 星历消息并按 PRN 存入星历数组
@param data 指向完整消息起始位置的字节缓冲区
@param geph GPS 星历记录数组，按 PRN-1 作为下标存放
@return 无
@note 各字段按消息体固定偏移逐字节读取，PRN 超出 1~MAXGPSNUM 时直接返回
*/
void ParseGpsEphemeris(unsigned char* data, EphemerisRecord geph[])
{
    unsigned char* p = data + 28;  // 跳过28字节头

    const int prn = (int)ReadAsUint32(p);  // PRN 在偏移 H+0
    if (prn < 1 || prn > MAXGPSNUM) return;

    EphemerisRecord* eph = geph + (prn - 1);

    eph->PRN = prn;
    eph->System = GPS;
    eph->SVHealth = ReadAsUint32(p + 12); // health

    // 时间与星历标识
    eph->TOE.Week = eph->TOC.Week = ReadAsUint32(p + 24); // 周编号
    eph->TOE.SecOfWeek = ReadAsDouble(p + 32);             // toe 秒
    eph->TOC.SecOfWeek = ReadAsDouble(p + 164);            // tow（subframe1 时间，先作为 TOC 秒用）
    eph->IODE = ReadAsUint32(p + 16);
    eph->IODC = ReadAsUint32(p + 20);

    // 钟差参数
    eph->ClkBias = ReadAsDouble(p + 180);       // a0
    eph->ClkDrift = ReadAsDouble(p + 188);      // a1
    eph->ClkDriftRate = ReadAsDouble(p + 196);  // a2
    eph->TGD1 = ReadAsDouble(p + 172);          // 组延迟 TGD

    // 平均角速度改正与轨道平面参数
    eph->DeltaN = ReadAsDouble(p + 48);
    eph->SqrtA = sqrt(ReadAsDouble(p + 40));
    eph->e = ReadAsDouble(p + 64);
    eph->M0 = ReadAsDouble(p + 56);
    eph->omega = ReadAsDouble(p + 72);

    // 升交点参数与轨道倾角
    eph->OMEGA = ReadAsDouble(p + 144);
    eph->OMEGADot = ReadAsDouble(p + 152);
    eph->i0 = ReadAsDouble(p + 128);
    eph->iDot = ReadAsDouble(p + 136);

    // 轨道摄动改正项
    eph->Cuc = ReadAsDouble(p + 80);
    eph->Cus = ReadAsDouble(p + 88);
    eph->Crc = ReadAsDouble(p + 96);
    eph->Crs = ReadAsDouble(p + 104);
    eph->Cic = ReadAsDouble(p + 112);
    eph->Cis = ReadAsDouble(p + 120);
    // Anti-spoofing、URA、N 暂不处理
}


/*!
@brief 解析 NovAtel OEM7 BDS 星历消息并按 PRN 存入星历记录数组
@param data 指向完整消息起始位置的字节缓冲区
@param beph BDS 星历记录数组，按 PRN-1 作为下标存放
@return 无
@note 含 TOE/TOC、钟差参数、TGD1 与轨道六根数，PRN 越界时直接返回
*/
void ParseBdsEphemeris(unsigned char* data, EphemerisRecord beph[])
{
    unsigned char* p = data + 28;  // 跳过同步头+消息头

    const int prn = (int)ReadAsUint32(p);
    if (prn < 1 || prn > MAXBDSNUM) return;

    EphemerisRecord* eph = beph + (prn - 1);

    eph->PRN = prn;
    eph->System = BDS;
    eph->SVHealth = ReadAsUint32(p + 16);

    // 时间与龄期
    eph->TOE.Week = eph->TOC.Week = ReadAsUint32(p + 4);
    eph->TOE.SecOfWeek = ReadAsUint32(p + 72);
    eph->TOC.SecOfWeek = ReadAsUint32(p + 40);
    eph->IODC = ReadAsUint32(p + 36);   // AODC
    eph->IODE = ReadAsUint32(p + 68);   // AODE

    // 钟差参数与 TGD
    eph->ClkBias = ReadAsDouble(p + 44);
    eph->ClkDrift = ReadAsDouble(p + 52);
    eph->ClkDriftRate = ReadAsDouble(p + 60);
    eph->TGD1 = ReadAsDouble(p + 20);

    // 轨道参数
    eph->SqrtA = ReadAsDouble(p + 76);
    eph->e = ReadAsDouble(p + 84);
    eph->omega = ReadAsDouble(p + 92);   // 近地点幅角
    eph->DeltaN = ReadAsDouble(p + 100);
    eph->M0 = ReadAsDouble(p + 108);
    eph->OMEGA = ReadAsDouble(p + 116);  // 升交点赤经
    eph->OMEGADot = ReadAsDouble(p + 124);
    eph->i0 = ReadAsDouble(p + 132);
    eph->iDot = ReadAsDouble(p + 140);

    // 轨道摄动改正
    eph->Cuc = ReadAsDouble(p + 148);
    eph->Cus = ReadAsDouble(p + 156);
    eph->Crc = ReadAsDouble(p + 164);
    eph->Crs = ReadAsDouble(p + 172);
    eph->Cic = ReadAsDouble(p + 180);
    eph->Cis = ReadAsDouble(p + 188);
}


/*!
@brief 解析 PSRPOS 消息中的经纬度高程、定位状态与精度信息
@param data 指向完整消息起始位置的字节缓冲区
@param pos 输出的单点定位结果结构体指针
@return 无
@note 经纬度高程单位分别为度、度、米
*/
void ParsePsrposMessage(unsigned char* data, PoseReport* pos)
{
    unsigned char* p = data + 28;   // 跳过28字节消息头

    // 大地坐标：纬度、经度、高程（高程 = 椭球高 + 大地水准面差距）
    pos->undulation = ReadAsFloat(p + 32);
    pos->Pos[1] = ReadAsDouble(p + 8);
    pos->Pos[0] = ReadAsDouble(p + 16);
    pos->Pos[2] = ReadAsDouble(p + 24) + pos->undulation;

    // 各方向标准差与三维标准差
    pos->lat_sigma = ReadAsFloat(p + 40);
    pos->lon_sigma = ReadAsFloat(p + 44);
    pos->hgt_sigma = ReadAsFloat(p + 48);
    pos->SigmaPos = sqrt(pos->lat_sigma * pos->lat_sigma
        + pos->lon_sigma * pos->lon_sigma
        + pos->hgt_sigma * pos->hgt_sigma);
}


/*!
@brief 在字节流中搜索 OEM7 同步头 AA 44 12
@param buff 字节缓冲区首地址
@param from 搜索起始下标
@param len 缓冲区有效字节数
@return 同步头起始下标，未找到时返回 len
@note 只在剩余长度足以容纳 3 字节同步头时才继续匹配
*/
static int FindSyncHeader(const unsigned char* buff, int from, int len)
{
    int i = from;

    while (i + 3 <= len)
    {
        if (buff[i] == 0xAA && buff[i + 1] == 0x44 && buff[i + 2] == 0x12)
        {
            return i;
        }
        i++;
    }
    return len;
}

/*!
@brief 按消息类型分发到对应的解码函数
@param msgId 消息标识
@param frame 指向整条消息起始位置（含 28 字节消息头）
@param obs 观测值输出
@param geph GPS 星历数组
@param beph BDS 星历数组
@param pos 单点定位结果输出
@return 1=本条为观测数据消息，0=其它类型消息
@note 未识别的消息不做处理，直接交给调用方跳过
*/
static int DispatchMessage(int msgId, unsigned char* frame, EpochObservation* obs,
    EphemerisRecord geph[], EphemerisRecord beph[], PoseReport* pos)
{
    switch (msgId)
    {
    case 43:                       // 观测数据消息
        ParseRangeMessage(frame, obs);
        return 1;

    case 7:                        // GPS 星历
        ParseGpsEphemeris(frame, geph);
        return 0;

    case 1696:                     // BDS 星历
        ParseBdsEphemeris(frame, beph);
        return 0;

    case 42:                       // 位置解
        ParsePsrposMessage(frame, pos);
        return 0;

    default:                       // 未识别的消息
        return 0;
    }
}

/*!
@brief 在字节流中搜索 OEM7 消息头，校验后分发到各解码函数
@param buff 待解析的字节缓冲区，解析后将剩余数据移回缓冲区首部
@param len 输入时为缓冲区字节数，输出时为未处理的剩余字节数
@param obs 输出的观测值结构体指针
@param geph 输出的 GPS 星历数组
@param beph 输出的 BDS 星历数组
@param pos 输出的单点定位结果结构体指针
@return 固定返回 0，表示本次解析流程结束
@note 同步头为 AA 44 12，CRC 校验失败则跳过 3 字节继续搜索；codemode=1 时处理一条观测数据即返回
*/
int ParseOem7Stream(unsigned char buff[], int& len, EpochObservation* obs, EphemerisRecord geph[], EphemerisRecord beph[], PoseReport* pos)
{
    int i = 0;          // 当前帧起始下标
    int Status = 0;     // 本次是否已解出观测数据
    int MsgLen = 0;     // 消息体长度
    int MsgID = 0;      // 消息标识

    while (i < len)
    {
        // 1、查找同步字符 AA 44 12
        i = FindSyncHeader(buff, i, len);
        if (i >= len) break;                 // 剩余字节不足 3 个

        // 2、检查是否有足够的字节读取完整的消息头（28字节）
        if (i + 28 > len) break;

        // 从消息头中解码消息长度 MsgLen 和消息类型 MsgID
        MsgID = ReadAsUint16(buff + i + 4);
        MsgLen = ReadAsUint16(buff + i + 8);

        // 3、检查整条消息是否完整（头28字节 + MsgLen 字节 + CRC 4字节）
        if (i + 28 + MsgLen + 4 > len) break;

        // 4、CRC检验：不通过则跳过当前同步头继续搜索
        const unsigned int crcRead = ReadAsUint32(buff + i + 28 + MsgLen);
        if (ComputeCrc32(buff + i, 28 + MsgLen) != crcRead)
        {
            i += 3;
            continue;
        }

        // 5、根据消息 ID 解码
        Status = DispatchMessage(MsgID, buff + i, obs, geph, beph, pos);

        // 6、跳过整条消息
        i += 28 + MsgLen + 4;

        // mode=1 时一条一条观测数据处理；mode=0 时持续处理所有数据
        if (Status == 1 && codemode == 1) break;
    }

    // 7、把未处理的剩余数据移到缓冲区首部
    if (len > i)
    {
        memmove(buff, buff + i, len - i);
    }
    len -= i;
    return 0;
}
