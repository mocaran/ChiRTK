#pragma once
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cfloat>
#include <vector>
#include<stdio.h>
#include<windows.h>
#include <math.h>
#include <memory.h>
#include <algorithm> 

using namespace Eigen;
using std::vector;
using namespace std;

/*-----------------------------------------------
    常数定义
------------------------------------------------*/

/* 导航卫星系统定义 */
enum GnssSystem { UnknownSystem = 0, GPS, BDS, GLONASS, GALILEO, QZSS };    // 默认下为位置系统UnknownSystem

/* 弧度-角度转换 */
#define PAI 3.14159265358979323846264338338279                  // 圆周率
#define Rad (PAI / 180.0)                                       // 角度转弧度
#define Deg (180.0 / PAI)                                       // 弧度转角度

/* WGS84坐标系常数 */
#define R_WGS84     6378137.0                                   // 地球半长轴 [米]
#define F_WGS84     (1.0 / 298.257223563)                       // 地球扁率
#define W_WGS84     7.2921151467e-5                             // 地球自转角速度 [弧度/秒]
#define GM_WGS84    3986005e+8                                  // 地心引力常数 [m^3/s^2]

/* CGCS2000坐标系常数 */
#define R_CGCS2000  6378137.0                                   // 地球半长轴 [米]
#define F_CGCS2000  (1.0 / 298.257222101)                       // 地球扁率
#define W_CGCS2000  7.292115e-5                                 // 地球自转角速度 [弧度/秒]
#define GM_CGCS2000 3986004.418e+8                              // 地心引力常数 [m^3/s^2]

/* 数据解码 */
#define C_Light 299792458.0                                     // 光速[米每秒]

#define MAXCHANNUM 36                                           // 接收机最多可同时处理的通道数量
#define MAXGPSNUM  32                                           // GPS系统最大卫星编号
#define MAXBDSNUM 63                                            // 北斗系统最大卫星编号
#define MAXRAWLEN 40960                                         // 原始数据缓冲区的最大长度（单位：字节）
#define POLYCRC32   0xEDB88320u                                 //校验码 

/* 卫星PVT解算 */
#define  FG1_GPS  1575.42E6                                     // L1信号频率
#define  FG2_GPS  1227.60E6                                     // L2信号频率
#define  FG1_BDS  1561.098E6                                    // B1信号的基准频率
#define  FG3_BDS  1268.520E6                                    // B3信号的基准频率

#define  WL1_GPS  (C_Light/FG1_GPS)                             // L1信号波长
#define  WL2_GPS  (C_Light/FG2_GPS)                             // L2信号波长
#define  WL1_BDS  (C_Light/FG1_BDS)                             // B1信号波长
#define  WL3_BDS  (C_Light/FG3_BDS)                             // B3信号波长

/* 误差改正 */
#define  Htrop      15000.0                                     // 对流层顶高度 m
#define  SeaLevel   0                                           // 海平面高度H0
#define  TempT0     15 + 273.16                                 // 温度T0
#define  AtmosPre   1013.25                                     // 气压p0
#define  RelHum     0.50                                        // 相对湿度RH0

#define  GFThres    0.05                                        // GF限差 cm
#define  MWThres    3.0                                         // MW限差 m

#define CFGITEMS 14                                             // 配置项的数量

extern int codemode;                                                // 模式选择，事后解算为1，实时解算为0

/*-----------------------------------------------
    数据结构
------------------------------------------------*/

/* GPS时间 (周, 周内秒) */
struct GpsTime
{
    unsigned short Week;
    double SecOfWeek;

    /*!
    @brief GpsTime 构造函数：把各成员初始化为 0 或无效值
    */
    GpsTime() : Week(0), SecOfWeek(0.0) {}
};

/* 大地坐标 (经度, 纬度, 高程) */
union GeodeticCoord
{
    struct {
        double latitude;    // 纬度
        double longitude;   // 经度
        double height;      // 高程
    };
    double blh[3];
};

/* 每颗卫星的观测数据  */
struct SatelliteObservation
{
    short Prn;          // 卫星编号
    GnssSystem System;     // 卫星所属系统

    double P[2];        // 伪距观测值（单位：米），最多两频
    double L[2];        // 载波相位观测值（单位：周），最多两频
    double D[2];        // 多普勒观测值（单位：Hz），最多两频

    double cn0[2];      // 载噪比 C/N0（单位：dB-Hz），最多两频
    double LockTime[2]; // 信号锁定时间（单位：秒），最多两频

    unsigned char half[2]; // 半周模糊度标志（通常为 0 或 1），最多两频

    bool Valid;         // 数据有效性标志（true为有效，false为无效）

    /*!
    @brief SatelliteObservation 构造函数：把各成员初始化为 0 或无效值
    */
    SatelliteObservation()
    {
        Prn = 0;
        System = UnknownSystem;
        for (int i = 0; i < 2; i++)
            P[i] = L[i] = D[i] = 0.0;
        cn0[0] = cn0[1] = 0.0;
        LockTime[0] = LockTime[1] = 0.0;
        half[0] = half[1] = 0.0;
        Valid = false;
    }
};

/*  卫星星历结构体   */
struct EphemerisRecord
{
    short PRN;         // 卫星编号
    GnssSystem System;    // 卫星系统

    GpsTime TOC;       // 钟差参考时间
    GpsTime TOE;       // 轨道参考时间

    double ClkBias;    // 卫星钟偏差（秒）
    double ClkDrift;   // 卫星钟漂移（秒/秒）
    double ClkDriftRate; // 卫星钟漂移率（秒/秒²）

    double IODE;       // 星历数据编号（轨道部分）
    double IODC;       // 星历数据编号（钟差部分）

    double SqrtA;      // 轨道长半轴的平方根（√a），单位：m^0.5
    double M0;         // 平近点角（弧度），参考时间的值
    double e;          // 轨道偏心率
    double OMEGA;      // 初始升交点经度（弧度）
    double i0;         // 轨道倾角（弧度）
    double omega;      // 升交距角（近地点幅角，弧度）

    double Crs;        // 正弦方向径向轨道改正项（米）
    double Cuc;        // 余弦方向轨道幅角改正项（弧度）
    double Cus;        // 正弦方向轨道幅角改正项（弧度）
    double Cic;        // 余弦方向升交点赤经改正项（弧度）
    double Cis;        // 正弦方向升交点赤经改正项（弧度）
    double Crc;        // 余弦方向径向轨道改正项（米）

    double DeltaN;     // 平均角速度改正值（弧度/秒）
    double OMEGADot;   // 升交点赤经变化率（弧度/秒）
    double iDot;       // 轨道倾角变化率（弧度/秒）

    int SVHealth;      // 卫星健康状态（0=健康，其他=不健康）

    double TGD1;       // 群延迟（秒）

    // ====== 构造函数 ======
    /*!
    @brief EphemerisRecord 构造函数：把各成员初始化为 0 或无效值
    */
    EphemerisRecord()
    {
        PRN = 0;
        System = UnknownSystem;  // 或 GnssSystem(0)，默认 未知

        TOC.Week = TOE.Week = 0;
        TOC.SecOfWeek = TOE.SecOfWeek = 0.0;

        ClkBias = 0.0;
        ClkDrift = 0.0;
        ClkDriftRate = 0.0;

        IODE = 0.0;
        IODC = 0.0;

        SqrtA = 0.0;
        M0 = 0.0;
        e = 0.0;
        OMEGA = 0.0;
        i0 = 0.0;
        omega = 0.0;

        Crs = 0.0;
        Cuc = 0.0;
        Cus = 0.0;
        Cic = 0.0;
        Cis = 0.0;
        Crc = 0.0;

        DeltaN = 0.0;
        OMEGADot = 0.0;
        iDot = 0.0;

        SVHealth = 0;

        TGD1 = 0.0;
    }

};

/* 每个历元的定位结果结构体定义 */
struct PoseReport
{
    double Pos[3];         // 接收机位置（经度、纬度、高程）[度, 度, 米]

    double SigmaPos;       // 三维位置标准差（m）

    float undulation;      // 大地水准面高度差（m）
    float lat_sigma;       // 纬度标准差（m）
    float lon_sigma;       // 经度标准差（m）
    float hgt_sigma;       // 高度标准差（m）

    /*!
    @brief PoseReport 构造函数：把各成员初始化为 0 或无效值
    */
    PoseReport()
    {
        Pos[0] = Pos[1] = Pos[2] = 0.0;

        SigmaPos = 0.0;

        undulation = 0.0f;
        lat_sigma = 0.0f;
        lon_sigma = 0.0f;
        hgt_sigma = 0.0f;
    }
};

/* 每颗卫星位置、速度和钟差等的中间计算结果 */
struct SatelliteState
{
    double SatPos[3];     // 卫星在ECEF坐标系下的三维位置 (X, Y, Z), 单位: 米
    double SatVel[3];     // 卫星在ECEF坐标系下的三维速度 (VX, VY, VZ), 单位: 米/秒
    double SatClkOft;     // 卫星钟差, 单位: 秒 (包含卫星钟偏移、漂移及相对论效应)
    double SatClkSft;     // 卫星钟速差, 单位: 无量纲 (用于精密定位中的钟速改正)
    double Elevation;     // 卫星仰角, 单位: 弧度 (范围: 0 ~ π/2)
    double Azimuth;       // 卫星方位角, 单位: 弧度 (范围: 0 ~ 2π, 从北方向顺时针计量)
    double TropCorr;      // 对流层延迟改正量, 单位: 米 (用于消除大气延迟对信号的影响)
    double Tgd1;          // 群延迟改正参数, 单位: 秒 (L1/L2频率间偏差)
    bool Valid;           // 数据有效性标志: 
    // false=没有可用星历或星历过期/不健康, 
    // true=卫星位置、钟差计算成功且可用

    /*!
    @brief SatelliteState 构造函数：把各成员初始化为 0 或无效值
    */
    SatelliteState()
    {
        // 初始化位置和速度为0
        SatPos[0] = SatPos[1] = SatPos[2] = 0.0;
        SatVel[0] = SatVel[1] = SatVel[2] = 0.0;

        // 初始仰角设为0
        Elevation = 0;

        // 初始化钟差、钟速差和改正参数为0
        SatClkOft = SatClkSft = 0.0;
        Tgd1 = TropCorr = 0.0;

        Azimuth = 0.0;
        // 默认标记为无效数据
        Valid = false;
    }
};

/* CombinationState组合观测 */
struct CombinationState
{
    short Prn;          // 卫星号
    GnssSystem Sys;        // 系统名称
    double MW, GF, PIF; // 组合观测值

    int n;              // 平滑计数

    /*!
    @brief CombinationState 构造函数：把各成员初始化为 0 或无效值
    */
    CombinationState()
    {
        Prn = n = 0;
        Sys = UnknownSystem;
        MW = GF = PIF = 0.0;
    }
};

/* 每个历元的观测数据定义 */
struct EpochObservation
{
    GpsTime    Time;                           // 当前历元的时间（GPS时间，包括周和周内秒）
    short      SatNum;                         // 当前历元观测到的卫星数量
    SatelliteObservation SatObs[MAXCHANNUM];             // 每颗卫星的观测数据（如伪距、载波等）
    SatelliteState     SatPVT[MAXCHANNUM];             // 每颗卫星的位置、速度、钟差计算结果
    CombinationState       ComObs[MAXCHANNUM];             // 每颗卫星的组合观测量（如MW组合、GF组合）

    /*!
    @brief EpochObservation 构造函数：把各成员初始化为 0 或无效值
    */
    EpochObservation()
    {
        SatNum = 0;
    }
};

/*  每颗卫星的单差观测数据定义  */
struct SingleDiffSatellite
{
    short    Prn;
    GnssSystem  System;
    short    Valid;         // 0代表无效
    double   dP[2], dL[2];  // 伪距和载波相位观测值单差
    short    nBase, nRover; // 存储单差观测值对应的基准和流动站的数值索引号
    unsigned char half[2];  // 半周模糊度标志（0表示有半周），两频

    /*!
    @brief SingleDiffSatellite 构造函数：把各成员初始化为 0 或无效值
    */
    SingleDiffSatellite()
    {
        Prn = nBase = nRover = 0;
        System = UnknownSystem;
        dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
        half[0] = 0; half[1] = 0;
        Valid = -1;
    }
};

/*  每个历元的单差观测数据定义  */
struct SingleDiffEpoch
{
    GpsTime    Time;                   // 当前单差观测对应的GPS时间（周内秒、周数等）
    short      SatNum;                 // 本历元参与单差的卫星数（有效配对数）
    SingleDiffSatellite   SdSatObs[MAXCHANNUM];   // 按卫星存放的单差观测（伪距/相位差、系统/PRN、有效标志等）
    CombinationState       SdCObs[MAXCHANNUM];     // 保存用于周跳检测的历史MW/GF组合数据（上一历元的参考）

    /*!
    @brief SingleDiffEpoch 构造函数：把各成员初始化为 0 或无效值
    */
    SingleDiffEpoch()
    {
        SatNum = 0;                    // 构造时默认单差卫星数为0
    }
};

/* 双差解算相关的数据结构定义 */
struct DoubleDiffEpoch
{
    int RefPrn[2], RefPos[2];                   // 参考星的PRN号与索引位置：[0]=GPS, [1]=BDS
    int Sats, DDSatNum[2];                      // 所有系统的双差卫星数，总数与各系统数量
    double FixedAmb[MAXCHANNUM * 4];            // 固定解的模糊度值（前AmbNum为最优，后为次优）
    double ResAmb[2], Ratio;                    // 模糊度残差与比值检验结果
    float  FixRMS;                              // 固定解的RMS误差（米）
    double dPos[3];                             // 计算得到的基线向量（ECEF）
    bool   bFixed;                              // 是否成功固定模糊度（true = 已固定）

    /*!
    @brief DoubleDiffEpoch 构造函数：把各成员初始化为 0 或无效值
    */
    DoubleDiffEpoch()
    {
        for (int i = 0; i < 2; i++) {
            DDSatNum[i] = 0;
            RefPrn[i] = RefPos[i] = -1;
        }
        Sats = 0;
        dPos[0] = dPos[1] = dPos[2] = 0.0;
        ResAmb[0] = ResAmb[1] = FixRMS = Ratio = 0.0;
        bFixed = false;
        for (int i = 0; i < MAXCHANNUM * 2; i++) {
            FixedAmb[2 * i + 0] = FixedAmb[2 * i + 1] = 0.0;
        }
    }
};

/* RTK浮点解数据定义 */
struct FloatSolution
{
    double dX[3];               // 基线向量（dx,dy,dz）Rov-Bas
    double N[MAXCHANNUM * 2];   // 双差模糊度
    double sigma;               // 验后单位权中误差
    double PDOP, HDOP, VDOP;    // 位置精度因子、水平精度因子、垂直精度因子
    double Qn[(MAXCHANNUM * 2) * (MAXCHANNUM * 2)]; // 模糊度的协因数矩阵 
    int stanum[2], totalsatnum; // 解算所用双差卫星数（不包括两个系统各自的基准星）-> 用于确定Qn的有效元素数量

    /*!
    @brief FloatSolution 构造函数：把各成员初始化为 0 或无效值
    */
    FloatSolution()
    {
        dX[0] = dX[1] = dX[2] = 0.0;
        PDOP = HDOP = VDOP = 0.0;
        sigma = 0.0;
        stanum[0] = stanum[1] = totalsatnum = 0;
        for (int i = 0; i < MAXCHANNUM * 2; i++)
        {
            N[i] = 0.0;
            for (int j = 0; j < MAXCHANNUM * 2; j++)
            {
                Qn[i * MAXCHANNUM * 2 + j] = 0.0;
            }
        }
    }
};

/* 每个历元的单点定位/测速结果及其精度指标 */
struct SinglePointResult
{
    GpsTime Time;                               // 当前历元时间
    double Position[3];                         // 解算出的ECEF坐标（单位：米）
    double Velocity[3];                         // 解算出的ECEF速度（单位：米/秒）
    double RcvClkOft[2];                        // 接收机钟差：[0]=GPS钟差, [1]=BDS钟差，单位：秒
    double PDOP, SigmaPos;                      // 定位精度因子、位置标准差（单位：米）
    short  GPSSatNum, BDSSatNum;                // 参与解算的GPS和BDS卫星数
    short  AllSatNum;                           // 所有参与观测的卫星数（总数）

    /*!
    @brief SinglePointResult 构造函数：把各成员初始化为 0 或无效值
    */
    SinglePointResult()
    {
        for (int i = 0; i < 3; i++)Position[i] = Velocity[i] = 0.0;
        RcvClkOft[0] = RcvClkOft[1] = 0.0;
        PDOP = SigmaPos = 999.9;
        GPSSatNum = BDSSatNum = AllSatNum = 0;
    }
};

/* RTK定位的数据结构定义 */
struct SolverData {
    EpochObservation   BasEpk_0;                            // 基准站的上一个观测数据（EpochObservation 结构体）
    EpochObservation   RovEpk_0;                            // 流动站的上一个观测数据（EpochObservation 结构体）
    EpochObservation   BasEpk_cur;                          // 基准站的当前观测数据（EpochObservation 结构体）
    EpochObservation   RovEpk_cur;                          // 流动站的当前观测数据（EpochObservation 结构体）
    SingleDiffEpoch SdObs;                               // 单差观测数据（由基准站与流动站的原始观测生成）
    DoubleDiffEpoch     DDObs;                               // 双差解算结果（包含模糊度固定信息、基线向量等）
    EphemerisRecord  GpsEph[MAXGPSNUM];                   // GPS系统的广播星历（每颗卫星一个 EphemerisRecord 结构）
    EphemerisRecord  BdsEph[MAXBDSNUM];                   // BDS系统的广播星历（每颗卫星一个 EphemerisRecord 结构）
    PoseReport BaseBest;                                // 基站的参考位置
    PoseReport RoverBest;                               // 流动站的参考位置
};

/* RTK定位后卡尔曼滤波结果 */
struct KalmanState
{
    GpsTime Time;
    double X[3 + MAXCHANNUM * 2], P[(3 + MAXCHANNUM * 2) * (3 + MAXCHANNUM * 2)];
    int Index[MAXCHANNUM], nSats;
    DoubleDiffEpoch DDObs, CurDDObs;           // 上一个历元和当前历元的双差观测值信息
    SingleDiffEpoch SDObs;                 // 上一个历元的单差观测值
    double X0[3 + MAXCHANNUM * 2], P0[(3 + MAXCHANNUM * 2) * (3 + MAXCHANNUM * 2)];  // 状态备份
    double ratio;                     // ratio值
    bool IsInit;                      // 滤波是否初始化
    bool IsFixed;                     // 滤波是否进行了LAMBDA固

    /*!
    @brief KalmanState 构造函数：把各成员初始化为 0 或无效值
    */
    KalmanState() {
        ratio = 0.0;
        IsInit = false;
        IsFixed = false;
        nSats = 0;
        for (int i = 0; i < MAXCHANNUM; i++) Index[i] = -1;
        for (int i = 0; i < 3 + MAXCHANNUM * 2; i++) {
            X[i] = X0[i] = 0.0;
            for (int j = 0; j < 3 + MAXCHANNUM * 2; j++) P[i * (3 + MAXCHANNUM * 2) + j] = P0[i * (3 + MAXCHANNUM * 2) + j] = 0.0;
        }
    }
};

/* 流动站配置信息 */
struct RoverConfig   // 配置信息
{
    int  IsFileData;                     //  1=文件, 0=实时,
    int RTKProcMode;                     // 1=EKF, 2=LSQ
    int    RovPort, RovBaud;             // 流动站COM端口设置
    char   BasNetIP[20], RovNetIP[20];   // 基站、流动站IP地址
    int BasNetPort, RovNetPort;          // 基站、流动站端口
    double CodeNoise, CPNoise;           // 伪距噪声
    double ElevThreshold;                // 高度角阈值
    double RatioThres;                   // Ratio检验阈值

    char  BasObsDatFile[256], RovObsDatFile[256];       //  观测数据的文件名
    char  FloatFile[256], FixedFile[256], SppFile[256], EkfFile[256]; //  结果数据文件名

    /*!
    @brief RoverConfig 构造函数：把各成员初始化为 0 或无效值
    */
    RoverConfig()
    {
        IsFileData = RTKProcMode = -1;
        RovPort = 0;  RovBaud = 0; BasNetPort = 0; RovNetPort = 0;
        CodeNoise = 0;
        CPNoise = 0;
        ElevThreshold = 0.0;
        RatioThres = 0.0;
    }
};

/*-----------------------------------------------
    坐标转换与误差计算 函数声明
------------------------------------------------*/
/*!
@brief 大地坐标转地心直角坐标，按给定椭球参数计算
@param geodetic 大地坐标[3]：[0]=经度、[1]=纬度（单位度）、[2]=高程（单位米）
@param ecef 输出地心直角坐标[3]，单位米
@param R 地球椭球长半轴，单位米
@param F 地球椭球扁率，无量纲
@return 无
@note 内部先把角度转弧度，再按卯酉圈曲率半径逐分量计算坐标
*/
void GeodeticToEcef(const double geodetic[3], double ecef[3], const double R, const double F);                      /* 大地坐标 -> 笛卡尔坐标 */
/*!
@brief 地心直角坐标转大地坐标，闭式求解纬度与高程
@brief 地心直角坐标转大地坐标，迭代求解纬度与高程
@param ecef 输入地心直角坐标[3]，单位米
@param geodetic 输出大地坐标[3]：[0]=纬度、[1]=经度（单位度）、[2]=高程（单位米）
@param F 地球椭球扁率，无量纲
@return 无
@note 对纬度做不动点迭代，初值由椭球改正项给出，收敛阈值1e-12，最多迭代1000次；不收敛时向标准错误输出警告
*/
void EcefToGeodetic(const double ecef[3], double geodetic[3], const double R, const double F);                      /* 笛卡尔坐标 -> 大地坐标 */
/*!
@brief 生成由地心直角坐标差转站心坐标的旋转矩阵
@param blh 测站大地坐标，取blh[0]为纬度、blh[1]为经度（单位度）
@param Rmat 输出3×3旋转矩阵，行依次对应东、北、天方向
@return 无
@note 直接左乘地心直角坐标差向量即可得到东、北、天分量
*/
void EnuRotationFromGeodetic(const GeodeticCoord* blh, Eigen::Matrix3d& Rmat);                                          /* 生成大地坐标 -> ENU 的旋转矩阵 */
/*!
@brief 求两点地心直角坐标差在站心坐标系下的东、北、天分量
@param Xs 目标点地心直角坐标[3]，单位米
@param Xr 测站地心直角坐标[3]，单位米
@param blh 测站大地坐标，用于生成站心旋转矩阵
@param dENU 输出坐标差[3]：东、北、天分量，单位米
@return 无
@note 内部先按Xs减Xr求差，再左乘EnuRotationFromGeodetic生成的旋转矩阵
*/
void EcefOffsetToEnu(const double Xs[], const double Xr[], const GeodeticCoord* blh, double dENU[]);              /* 计算 ENU 坐标差（误差） */
/*!
@brief 计算卫星相对测站的高度角与方位角
@param Xr 测站地心直角坐标[3]，单位米
@param Xs 卫星地心直角坐标[3]，单位米
@param blh 测站大地坐标，用于生成站心旋转矩阵
@param Elev 输出高度角，单位度，范围0~90
@param Azim 输出方位角，单位度，由北方向顺时针计量，范围0~360
@return 无
@note 先求站心东、北、天坐标，再用atan2求角度，方位角为负时加360度
*/
void SatelliteElevationAzimuth(const double Xr[], const double Xs[], const GeodeticCoord* blh, double* Elev, double* Azim); /* 计算卫星高度角与方位角 */

/*-----------------------------------------------
    数据解码 函数声明
------------------------------------------------*/
/*!
@brief 从原始字节流中搜索并逐条解码NovAtel OEM7格式消息
@param buff 输入输出缓冲区，函数结束时把未处理的剩余字节移到缓冲区头部
@param len 输入缓冲区有效字节数，输出剩余未处理字节数
@param obs 输出当前历元观测数据，收到观测消息时填充
@param geph GPS星历数组，收到GPS星历消息时按卫星号写入
@param beph BDS星历数组，收到BDS星历消息时按卫星号写入
@param pos 输出单点定位结果，收到PSRPOS消息时填充
@return 固定返回0，表示本次处理结束，不是错误码
@note 按同步头AA 44 12定位消息并做CRC32校验；codemode为1时解出一条观测消息即返回，否则处理完整个缓冲区
*/
int ParseOem7Stream(unsigned char buff[], int& len, EpochObservation* obs, EphemerisRecord geph[], EphemerisRecord beph[], PoseReport* pos);   /* 解码主函数 */
/*!
@brief 计算指定长度字节流的CRC32校验码
@param buff 待校验数据的首地址
@param len 参与校验的字节数
@return 32位CRC校验值，用于与消息尾部的校验字段比较
@note 采用多项式0xEDB88320，按位右移逐位计算，不查表
*/
unsigned int ComputeCrc32(const unsigned char* buff, int len);             /* CRC32校验码 */
/*!
@brief 解码OEM7二进制观测数据消息，填充当前历元的卫星观测值
@param msg 指向一条完整消息的起始地址，含28字节消息头
@param obs 输出历元观测数据，写入观测时刻、各卫星伪距相位多普勒与卫星数
@return 无
@note 只保留GPS的L1/L2与BDS的B1I/B3I信号；伪距或相位未锁定时对应观测值置0，载波相位单位周
*/
void ParseRangeMessage(const unsigned char* msg, EpochObservation* obs);          /* 解码二进制观测数据消息 */
/*!
@brief 解码GPS广播星历消息并按卫星号存入星历数组
@param msg 指向一条完整GPS星历消息的起始地址
@param geph GPS星历数组，按卫星号减1的下标写入该卫星的星历参数
@return 无
@note 卫星号超出1~32范围时直接返回；消息中的√a开方后存入SqrtA
*/
void ParseGpsEphemeris(const unsigned char* msg, EphemerisRecord geph[]);    /* 解码GPS观测值数据 */
/*!
@brief 解码BDS广播星历消息并按卫星号存入星历数组
@param msg 指向一条完整BDS星历消息的起始地址
@param beph BDS星历数组，按卫星号减1的下标写入该卫星的星历参数
@return 无
@note 卫星号超出1~63范围时直接返回；TGD1解码
*/
void ParseBdsEphemeris(const unsigned char* msg, EphemerisRecord beph[]);    /* 解码BDS观测值数据 */
/*!
@brief 解码NovAtel单点定位解消息，填充定位结果结构体
@param msg 指向一条完整PSRPOS消息的起始地址
@param pos 输出定位结果，写入解状态、纬度经度高程、卫星数与各方向标准差
@return 无
@note 三维位置标准差由纬、经、高三个方向标准差的平方和开方得到，速度标准差与PDOP置0
*/
void ParsePsrposMessage(const unsigned char* msg, PoseReport* pos);          /* 解码单点定位解 */

/*-----------------------------------------------
    定位解算 函数声明
------------------------------------------------*/
/*!
@brief 判断星历是否可用并计算不含相对论改正的卫星钟差
@param Prn 卫星编号，作为星历数组的下标使用
@param Sys 卫星系统，仅支持GPS与BDS
@param t 当前观测时刻的GPS时间
@param GPSEph GPS星历数组
@param BDSEph BDS星历数组
@param Mid 卫星中间结果，输出SatClkOft卫星钟差，单位秒
@return true=星历有效且钟差已算出 false=系统不支持、星历过期或卫星不健康
@note GPS星历有效期7500秒、BDS为3900秒；计算BDS钟差前先把时间减去1356周与14秒
*/
bool SatelliteClockCorrection(const int Prn, const GnssSystem Sys, const GpsTime* t, EphemerisRecord* GPSEph, EphemerisRecord* BDSEph, SatelliteState* Mid);  /* 星历是否过期的判断    不加相对论改正的卫星钟差 */
/*!
@brief 由GPS广播星历计算卫星位置、速度、钟差与钟速
@param Prn 卫星编号
@param t 信号发射时刻的GPS时间
@param Eph 指向该卫星的星历记录
@param Mid 输出卫星位置速度、钟差钟速、TGD与有效性标志
@return 0=计算成功 1=星历过期或卫星不健康，此时只把Valid置为false
@note 钟差含相对论效应改正，位置速度在WGS84地固系下给出，单位米与米每秒
*/
int ComputeGpsSatelliteState(const int Prn, const GpsTime* t, const EphemerisRecord* Eph, SatelliteState* Mid);                                      /* 计算GPS卫星位置、速度、钟差 */
/*!
@brief 由BDS广播星历计算卫星位置、速度、钟差与钟速
@param Prn 卫星编号
@param t 信号发射时刻的GPS时间，函数内部换算为北斗时
@param Eph 指向该卫星的星历记录
@param Mid 输出卫星位置速度、钟差钟速、TGD与有效性标志
@return 0=计算成功 1=星历过期或卫星不健康，此时只把Valid置为false
@note 轨道倾角小于30度时按GEO卫星处理，需做绕X轴-5度与绕Z轴的坐标旋转
*/
int ComputeBdsSatelliteState(const int Prn, const GpsTime* t, const EphemerisRecord* Eph, SatelliteState* Mid);                                      /* 计算BDS卫星位置、速度、钟差 */

/*-----------------------------------------------
    误差改正 函数声明
------------------------------------------------*/
/*!
@brief 用HopfieldTroposphereDelay模型计算对流层延迟改正量
@param H 测站高程，单位米
@param Elev 卫星高度角，单位度
@return 对流层延迟改正量，单位米；测站高程超过对流层顶15000米时返回0
@note 干、湿分量分别按各自等效高度计算，再经含高度角的映射函数投影到传播路径
*/
double HopfieldTroposphereDelay(const double H, const double Elev);                                                             /* 对流层延迟改正 */
/*!
@brief 用GF与MW组合的历元间变化探测伪距相位粗差
@param Obs 历元观测数据，读取双频伪距相位，就地更新Valid与组合观测缓存ComObs
@return 无
@note GF限差0.05米、MW限差3米，任一超限即判为粗差；首个历元无历史组合值时卫星保持无效
*/
void RejectObservationOutliers(EpochObservation* Obs);                                                                              /* 粗差探测 */

/*-----------------------------------------------
    SolveSinglePointPosition & SolveSinglePointVelocity 函数声明
------------------------------------------------*/
/*!
@brief 迭代求信号发射时刻并计算该时刻各卫星的位置速度与钟差
@param Epk 历元观测数据，提供观测值与接收机近似位置，输出各卫星的SatPVT
@param GPSEph GPS星历数组
@param BDSEph BDS星历数组
@param UserPos 接收机近似坐标[3]，单位米，用于计算传播时间与地球自转改正
@return 无
@note 钟差按1e-12秒阈值迭代，随后做地球自转改正并计算高度角、方位角与对流层延迟
*/
void SatelliteStateAtSignalTime(EpochObservation* Epk, EphemerisRecord* GPSEph, EphemerisRecord* BDSEph, double UserPos[3]);    /*	信号发射时刻卫星位置计算 */
/*!
@brief 用双频无电离层组合伪距做单点定位迭代最小二乘解算
@param Epoch 当前历元观测数据，提供观测值、卫星位置速度与组合观测
@param Raw 原始数据结构体，提供GPS与BDS星历数组
@param Result 定位结果，输入上一历元位置作初值，输出坐标、钟差、DOP与精度
@return true=迭代收敛并解出坐标 false=有效卫星数少于待估参数个数
@note 待估参数为三维坐标加各系统接收机钟差，BDS另加Tgd1改正；按系统标记紧凑组织设计阵，最多迭代16次，位置改正量小于1e-7米判为收敛
*/
bool SolveSinglePointPosition(EpochObservation* Epoch, SolverData* Raw, SinglePointResult* Result);                                                   /*  单点定位SolveSinglePointPosition解算 */
/*!
@brief 用多普勒观测做单点测速最小二乘解算
@param Epoch 当前历元观测数据，提供多普勒观测值与卫星位置速度
@param Result 测速结果，输入已解出的Position，输出Velocity与速度标准差
@return 无
@note 有效卫星少于4颗或法方程矩阵不可逆时不更新结果并直接返回
*/
void SolveSinglePointVelocity(EpochObservation* Epoch, SinglePointResult* Result);                                                                /*  单点测速SolveSinglePointVelocity解算 */

/*-----------------------------------------------
    Socket网络
------------------------------------------------*/
/*!
@brief 以TCP方式连接指定IP与端口的GNSS数据流服务
@param sock 输出套接字句柄，失败时置为INVALID_SOCKET
@param IP 服务器IP地址字符串，如192.168.1.100
@param Port 服务器端口号
@return true=连接成功 false=初始化、建套接字、连接或设置接收超时失败
@note 连接成功后把接收超时设为5000毫秒，失败路径会关闭套接字并清理Winsock
*/
bool ConnectStreamSocket(SOCKET& sock, const char IP[], const unsigned short Port);  /*  打开串口通信  */
/*!
@brief 关闭已打开的套接字并释放Winsock资源
@param sock 待关闭的套接字句柄
@return 无
@note 内部调用closesocket与WSACleanup，调用后不应再使用该句柄
*/
void DisconnectStreamSocket(SOCKET& sock); /*  关闭串口通信  */

/*-----------------------------------------------
    内置配置解析
------------------------------------------------*/
/*!
@brief 逐行解析编译进程序的内置配置文本并解析14个配置项到结构体
@param ConfigName 内置配置的名称（如 LSQ_FILE），仅用于日志与错误提示
@param ANinfo 输出的流动站配置信息，写入运行模式、端口、文件名与各项阈值（配置文本来自 builtin_config.h）
@return true=全部配置项解析成功 false=文本为空或任一配置项格式非法
@note 感叹号!之后的内容按注释截断；文件模式与实时模式分别解析不同的配置项
*/
bool ParseRoverConfig(const char ConfigName[], const char ConfigText[], RoverConfig* ANinfo);   /* 解析内置配置 */

/*-----------------------------------------------
    RTK函数
------------------------------------------------*/
/*!
@brief 从基站与流动站观测文件中各读一个历元并按时标同步
@param FBase 基站观测文件指针
@param FRover 流动站观测文件指针
@param Raw 原始数据，输出同步后的基站历元BasEpk与流动站历元RovEpk
@return 1=同步成功 0=时刻不匹配或流动站滞后 -1=文件读取失败或已到文件尾
@note 时间差限差0.01秒；同步前把当前历元备份到BasEpk_0与RovEpk_0用于周跳探测
*/
int SyncEpochFromFiles(FILE* FBase, FILE* FRover, SolverData* Raw);                                        /* 基站与流动站数据时间同步（文件）*/
/*!
@brief 从基站与流动站网络端口中各收一个历元并按时标同步
@param BaseCom 基站通信套接字
@param RoverCom 流动站通信套接字
@param Raw 原始数据，输出同步后的基站历元BasEpk与流动站历元RovEpk
@return 1=同步成功 0=未同步或接收超时 -2=网络关闭或读取失败
@note 时间差限差5秒；每次调用先休眠约1秒等待新数据到达
*/
int SyncEpochFromSockets(SOCKET& BaseCom, SOCKET& RoverCom, SolverData* Raw);                                /* 基站与流动站数据时间同步（实时）*/
/*!
@brief 按连续跟踪时长标记各卫星观测是否可用于解算
@param EpkLast 上一历元观测数据，提供各卫星的失锁与跟踪时长信息
@param EpkCur 当前历元观测数据，就地更新各卫星的Valid标志
@return 无
@note 上一历元缺失或失锁时要求本历元LockTime不小于6秒，双频观测都通过才标记有效
*/
void MarkObservationsUsable(const EpochObservation* EpkLast, EpochObservation* EpkCur);                                      /* 观测数据有效性的标记 */
/*!
@brief 对基站与流动站同一卫星的观测做站间单差
@param EpkR 流动站历元观测数据
@param EpkB 基站历元观测数据
@param SDObs 输出站间单差观测，含伪距、载波相位与半周标志
@return 无
@note 仅两站同一卫星都有效时才做单差，缺失频率置0；半周标志要求两站均正常才置1
*/
void FormInterStationDifference(const EpochObservation* EpkR, const EpochObservation* EpkB, SingleDiffEpoch* SDObs);                   /* 基站与流动站 站间单差 */
/*!
@brief 用GF与MW组合探测站间单差观测中的周跳
@param Obs 站间单差观测，就地更新Valid标志与组合观测平滑缓存
@return 无
@note 只用GPS与BDS；GF限差0.05米、MW限差3米，超限即判为周跳
*/
void DetectCycleSlipEvents(SingleDiffEpoch* Obs);                                                          /* 周跳探测函数 */
/*!
@brief 按高度角与信噪比为各系统选取双差参考星
@param EpkA 流动站历元观测数据
@param EpkB 基站历元观测数据
@param SDObs 站间单差观测数据
@param DDObs 双差观测数据，输出GPS与BDS的参考星编号RefPrn
@return 无
@note 综合质量分需大于200且高度角大于40度，否则把RefPrn置为-1表示该系统无参考星
*/
void ChooseReferenceSatellite(const EpochObservation* EpkA, const EpochObservation* EpkB, SingleDiffEpoch* SDObs, DoubleDiffEpoch* DDObs);   /* 选取参考卫星 */
/* 最小二乘法 */
/*!
@brief 用双差伪距与载波相位做浮点解迭代最小二乘
@param Rawdata 原始数据，提供同步历元、单差双差观测与配置参数
@param Base 基站单点定位结果
@param Rover 流动站单点定位结果
@param Fres 输出浮点解的位置改正dX、模糊度N、协因数Qn、DOP与单位权中误差
@return true=解算收敛 false=双差卫星少于5颗或迭代50次仍未收敛
@note 每颗卫星4个观测方程，相位权重取伪距权重的1000倍
*/
bool EstimateFloatSolution(SolverData* Rawdata, SinglePointResult* Base, SinglePointResult* Rover, FloatSolution* Fres);             /* 相对定位浮点解函数设计 */
/*!
@brief 从法方程逆阵中提取模糊度参数的协因数子块
@param satnum 参与解算的双差卫星数
@param N_inv 法方程系数阵的逆矩阵
@param Q 输出模糊度协因数矩阵，按行优先存放
@return 无
@note 前3个参数为位置改正提取时跳过，其后每颗卫星对应双频两个模糊度参数
*/
void ExtractAmbiguityCofactor(const int satnum, const MatrixXd& N_inv, double* Q);                                /* 提取QNN协因数矩阵 */
/*!
@brief 固定模糊度后只用双差相位解算位置改正
@param Rawdata 原始数据，提供双差观测与已固定的整周模糊度
@param Base 基站单点定位结果
@param Rover 流动站单点定位结果
@param Fres 输出固定解的位置改正与精度信息
@return true=解算收敛 false=迭代50次仍未收敛
@note 只估计3个位置改正，模糊度直接取Rawdata->DDObs.FixedAmb中的整数值
*/
bool EstimateFixedBaseline(SolverData* Rawdata, SinglePointResult* Base, SinglePointResult* Rover, FloatSolution* Fres);               /* 相对定位固定解 */
/* 卡尔曼滤波法 */
/*!
@brief 初始化卡尔曼滤波的状态向量、协方差与卫星索引
@param Rawdata 原始数据，提供双差观测与卫星编号
@param Base 基站单点定位结果
@param Rover 流动站单点定位结果
@param ekf 滤波结构体，输出状态X0、协方差P0、卫星索引Index与初始化标志
@return 无
@note 状态由3个坐标改正加双频双差模糊度组成，坐标初值方差100、模糊度500，参考星与无效卫星索引置-1
*/
void ResetKalmanFilter(SolverData* Rawdata, SinglePointResult* Base, SinglePointResult* Rover, KalmanState* ekf);                 /* 滤波初始化 */
/*!
@brief 卡尔曼滤波时间更新，把状态与协方差递推到当前历元
@param Rawdata 原始数据，提供当前历元的双差卫星与参考星信息
@param ekf 滤波结构体，就地更新状态X与协方差P
@return true=递推成功 false=状态维数与当前卫星不匹配导致越界
@note 位置按静态处理、过程噪声取1e-8；参考星变化时用新旧参考星之差传递状态，无对应卫星按方差100与500重新初始化
*/
bool PredictKalmanState(SolverData* Rawdata, KalmanState* ekf);                                                  /* 时间更新函数 */
/*!
@brief 求双差观测设计矩阵行，即卫星与参考星方向余弦之差
@param Rover 流动站坐标[3]，单位米
@param Sat 卫星坐标[3]，单位米
@param Ref 参考星坐标[3]，单位米
@return 卫星视线向量与参考星视线向量之差构成的3维向量，即设计矩阵H的前三列
@note 视线向量为单位向量，双差后两星公共的接收机钟差项被消去
*/
inline Vector3d DoubleDifferenceUnitVector(const double Rover[3], const double Sat[3], const double Ref[3]); /* 方向余弦计算 */
/*!
@brief 卡尔曼滤波量测更新，用双差观测修正状态并检核模糊度
@param Rawdata 原始数据，提供双差观测值、配置中的噪声与ratio阈值
@param ekf 滤波结构体，就地更新状态、协方差、ratio与固定标志
@return true=量测更新成功 false=有效观测不足或矩阵运算失败
@note 每颗卫星双频伪距加相位共4个观测，噪声取CodeNoise与CPNoise；ratio达到RatioThres后二次更新固定模糊度，否则置IsInit为false等待重新初始化
*/
bool CorrectKalmanState(SolverData* Rawdata, KalmanState* ekf);                                            /* 测量更新函数 */

/*-----------------------------------------------
    Lambda函数
------------------------------------------------*/
/*!
@brief 用LAMBDA方法对双差模糊度做整数最小二乘固定
@param n 待固定模糊度的维数
@param m 模糊度参数个数，通常与n相同
@param a 浮点模糊度向量，长度n
@param Q 模糊度协方差矩阵，按列优先存放，尺寸n×n
@param F 输出的固定解与检验量
@param s 输出的整数最小二乘残差平方和
@return 0=固定成功 -1=维数非法或分解、搜索失败
@note 流程为先做CholeskyLdDecomposition分解，再约化基向量并搜索，最后把结果逆变换回原模糊度空间
*/
int SearchIntegerLeastSquares(int n, int m, const double* a, const double* Q, double* F, double* s);                       /* SearchIntegerLeastSquares/mlambda整数最小二乘估计 */
/*!
@brief 用高斯-约当消元（列主元）求n阶方阵的逆矩阵
@param n 矩阵阶数
@param a 输入的n×n方阵，按行优先存放，函数内不修改
@param b 输出的逆矩阵，按行优先存放，尺寸与a相同
@return 1=求逆成功 0=阶数非法或主元素接近0导致矩阵不可逆
@note 逐列选取主元并记录行列交换，主元素小于DBL_EPSILON时判为不可逆；消元后按记录逆序恢复行列次序
*/
int InvertMatrix(int n, double a[], double b[]);                                                           /* 矩阵求逆 */
/*!
@brief 计算两个按行优先存放的矩阵乘积
@param m1 左矩阵行数
@param n1 左矩阵列数，同时作为右矩阵行数的约束
@param m2 右矩阵行数，必须与n1相等
@param n2 右矩阵列数
@param M1 左矩阵元素，按行优先存放，长度m1×n1
@param M2 右矩阵元素，按行优先存放，长度m2×n2
@param M3 输出乘积矩阵，按行优先存放，长度m1×n2
@return true=相乘成功 false=维数非正或n1与m2不相等
@note 采用i-j-k三重循环，内层沿k累加，直接读取两输入矩阵；输出矩阵不可与输入矩阵共用内存
*/
bool MultiplyMatrices(int m1, int n1, int m2, int n2, const double M1[], const double M2[], double M3[]); /* 矩阵的乘法 */
/*!
@brief 在约化后的基向量上搜索整数最小二乘解
@param n 模糊度维数
@param m 候选整数的个数
@param L 下三角变换矩阵，按列优先存放
@param D 对角阵元素，长度n
@param zs 变换后的浮点解向量，长度n
@param zn 输出的整数解向量，长度n
@param s 输出的残差平方和，最优值保存在s[0]
@return 0=搜索成功 -1=超过最大循环次数LOOPMAX仍未找到可行解
@note 采用深度优先收缩搜索，从最后一维逐层向内枚举候选整数
*/
int SearchIntegerCandidates(int n, int m, const double* L, const double* D, const double* zs, double* zn, double* s);    /* 改进的SearchIntegerLeastSquares搜索算法 */
/*!
@brief 对基向量做整数高斯约化以改善搜索效率
@param n 模糊度维数
@param L 待约化的下三角矩阵，按列优先存放，就地更新
@param D 对角阵元素，长度n，就地更新
@param Z 输出的整数变换矩阵，累积约化过程中的变换
@return 无
@note 交换相邻列时附加1E-6的容差判断，避免因舍入误差反复交换
*/
void LllBasisReduction(int n, double* L, double* D, double* Z);                                                 /* SearchIntegerLeastSquares约化 */
/*!
@brief 交换基向量中第j与第j+1列并同步更新变换矩阵
@param n 模糊度维数
@param L 下三角矩阵，就地交换相关的列
@param D 对角阵元素，就地交换对应元素
@param j 待交换列的起始下标
@param del 交换后对角元素的增量
@param Z 整数变换矩阵，同步交换对应的列
@return 无
@note 本函数是约化过程的内部步骤，交换后需重新计算L与D的相应元素
*/
void SwapBasisColumns(int n, double* L, double* D, int j, double del, double* Z);                                   /* 置换变换 */
/*!
@brief 用整数高斯变换消去第j列对第i列的影响
@param n 模糊度维数
@param L 下三角矩阵，就地做整数消元
@param Z 整数变换矩阵，同步做相同的整数变换
@param i 被消元的目标列下标
@param j 参与消元的源列下标
@return 无
@note 变换系数由mu=ROUND(L[i+j*n])取整得到，保证变换矩阵仍为整数矩阵
*/
void GaussIntegerTransform(int n, double* L, double* Z, int i, int j);                                                  /* 整数高斯变换 */

/*-----------------------------------------------
    结果打印函数
------------------------------------------------*/
/*!
@brief 把单点定位结果按历元逐行写入结果文件
@param rawdata 原始观测数据，取流动站当前历元的时间
@param Fres 单点定位结果，提供坐标、DOP、位置标准差与各系统卫星数
@param enu 定位误差的东、北、天分量[3]，单位米
@param OutputFile 已打开的结果文件指针
@return 无
@note 每行依次输出周、周内秒、ECEF坐标、ENU误差、PDOP、位置标准差与BDS/GPS卫星数
*/
void WriteSinglePointRecord(SolverData* rawdata, SinglePointResult* Fres, double enu[3], FILE* OutputFile);                             /* 输出单点定位结果 */
/*!
@brief 把相对定位浮点解结果按历元逐行写入结果文件
@param rawdata 原始观测数据，取流动站当前历元时间与双差卫星数
@param ecef 流动站坐标[3]，单位米
@param denu 相对基站位置误差的东、北、天分量[3]，单位米
@param Fres 浮点解结果，提供单位权中误差与各方向DOP
@param OutputFile 已打开的结果文件指针
@return 无
@note 每行依次输出时间、ECEF坐标、ENU误差、单位权中误差、PDOP/HDOP/VDOP与两系统双差卫星数
*/
void WriteFloatRecord(SolverData* rawdata, double ecef[3], double denu[3], FloatSolution* Fres, FILE* OutputFile);        /* 输出浮点解定位结果 */
/*!
@brief 把相对定位固定解结果按历元逐行写入结果文件
@param rawdata 原始观测数据，取流动站当前历元时间与双差卫星数
@param ecef 流动站坐标[3]，单位米
@param denu 相对基站位置误差的东、北、天分量[3]，单位米
@param OutputFile 已打开的结果文件指针
@return 无
@note 每行依次输出时间、ECEF坐标、ENU误差、ratio值与两系统双差卫星数
*/
void WriteFixedRecord(SolverData* rawdata, double ecef[3], double denu[3], FILE* OutputFile);                           /* 输出固定解定位结果 */
/*!
@brief 把卡尔曼滤波相对定位结果按历元逐行写入结果文件
@param rawdata 原始观测数据，取流动站当前历元的时间
@param ekf 滤波结构体，取状态X0前三个元素作为坐标、ratio作为检核值
@param denu 相对基站位置误差的东、北、天分量[3]，单位米
@param OutputFile 已打开的结果文件指针
@return 无
@note 每行依次输出时间、ECEF坐标、ENU误差、ratio值与两系统双差卫星数
*/
void WriteKalmanRecord(SolverData* rawdata, KalmanState* ekf, double denu[3], FILE* OutputFile);                               /* 输出卡尔曼定位结果 */
/*!
@brief 按RTKPLOT格式输出定位结果用于绘图分析
@param rawdata 原始观测数据，取流动站当前历元的时间
@param FixedBlh 最小二乘固定解的大地坐标[3]：纬度/经度/高程，单位度、度、米
@param ekf 滤波结构体，模式1时提供滤波坐标、ratio与固定标志
@param EkfBlh 卡尔曼滤波解的大地坐标[3]：纬度/经度/高程，单位度、度、米
@param OutputFile 已打开的结果文件指针
@param mode 输出模式：1=卡尔曼滤波解，2=最小二乘解
@return 无
@note 解状态1=固定解、2=浮点解，由模糊度是否固定决定；每行输出周、周内秒、BLH、ratio与解状态
*/
void WritePlotRecord(SolverData* rawdata, double FixedBlh[3], KalmanState* ekf, double EkfBlh[3], FILE* OutputFile, int mode); /* 输出定位结果用于RTKPLOT画图 */
