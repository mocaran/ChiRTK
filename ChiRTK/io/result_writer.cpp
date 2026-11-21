#include "gnss_core_types.h"


namespace
{
    /*!
    @brief 写出一行结果的公共前 8 列：周、周内秒、三维坐标、东北天偏差
    @param OutputFile 已打开的结果文件句柄
    @param Time 当前历元的 GPS 时间
    @param pos 三维坐标数组，单位米（或大地坐标下的度、度、米）
    @param denu 东、北、天方向的位置偏差数组，单位米
    @return 无
    @note 各结果文件的前 8 列格式统一在此维护，后续各函数只需追加各自的尾部列
    */
    void WriteEpochCoordinatePrefix(FILE* OutputFile, const GpsTime& Time, const double pos[3], const double denu[3])
    {
        fprintf(OutputFile, "%u %.f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f",
            Time.Week, Time.SecOfWeek,
            pos[0], pos[1], pos[2],
            denu[0], denu[1], denu[2]);
    }
}

/*!
@brief 将单点定位结果按固定列格式写入 SolveSinglePointPosition 结果文件
@param rawdata 原始观测数据指针，用于获取当前历元时间
@param Fres 单点定位结果结构体指针
@param enu 东、北、天方向的位置偏差数组，单位米
@param OutputFile 已打开的结果文件句柄
@return 无
*/
void WriteSinglePointRecord(SolverData* rawdata, SinglePointResult* Fres, double enu[3], FILE* OutputFile)
{
    WriteEpochCoordinatePrefix(OutputFile, rawdata->RovEpk_cur.Time, Fres->Position, enu);

    fprintf(OutputFile, " %0.15f %0.15f %d %d\n",
        Fres->PDOP, Fres->SigmaPos,
        Fres->BDSSatNum, Fres->GPSSatNum);
}

/*!
@brief 将浮点解定位结果按固定列格式写入浮点解结果文件
@param rawdata 原始观测数据指针，用于获取当前历元时间与差分卫星数
@param ecef 地心直角坐标数组，单位米
@param denu 东、北、天方向的位置偏差数组，单位米
@param Fres 浮点解结果结构体指针
@param OutputFile 已打开的结果文件句柄
@return 无
*/
void WriteFloatRecord(SolverData* rawdata, double ecef[3], double denu[3], FloatSolution* Fres, FILE* OutputFile)
{
    WriteEpochCoordinatePrefix(OutputFile, rawdata->RovEpk_cur.Time, ecef, denu);

    fprintf(OutputFile, " %0.8f %0.3f %0.3f %0.3f %d %d\n",
        Fres->sigma, Fres->PDOP, Fres->HDOP, Fres->VDOP,
        rawdata->DDObs.DDSatNum[0], rawdata->DDObs.DDSatNum[1]);
}

/*!
@brief 将固定解定位结果按固定列格式写入固定解结果文件
@param rawdata 原始观测数据指针，用于获取当前历元时间与 Ratio 值
@param ecef 地心直角坐标数组，单位米
@param denu 东、北、天方向的位置偏差数组，单位米
@param OutputFile 已打开的结果文件句柄
@return 无
*/
void WriteFixedRecord(SolverData* rawdata, double ecef[3], double denu[3], FILE* OutputFile)
{
    WriteEpochCoordinatePrefix(OutputFile, rawdata->RovEpk_cur.Time, ecef, denu);

    fprintf(OutputFile, " %0.3f %d %d\n",
        rawdata->DDObs.Ratio,
        rawdata->DDObs.DDSatNum[0], rawdata->DDObs.DDSatNum[1]);
}

/*!
@brief 将卡尔曼滤波定位结果按固定列格式写入 EKF 结果文件
@param rawdata 原始观测数据指针，用于获取当前历元时间
@param ekf 卡尔曼滤波状态结构体指针，提供位置与 Ratio 值
@param denu 东、北、天方向的位置偏差数组，单位米
@param OutputFile 已打开的结果文件句柄
@return 无
*/
void WriteKalmanRecord(SolverData* rawdata, KalmanState* ekf, double denu[3], FILE* OutputFile)
{
    WriteEpochCoordinatePrefix(OutputFile, rawdata->RovEpk_cur.Time, ekf->X0, denu);

    fprintf(OutputFile, " %0.3f %d %d\n",
        ekf->ratio, ekf->CurDDObs.DDSatNum[0], ekf->CurDDObs.DDSatNum[1]);
}

// mode=1是EKF，mode=2是LSQ
/*!
@brief 按处理模式输出用于绘制 RTK 轨迹图的历元定位结果
@param rawdata 原始观测数据指针，用于获取当前历元时间与固定标志
@param FixedBlh 最小二乘固定解的大地坐标数组，单位度、度、米
@param ekf 卡尔曼滤波状态结构体指针，mode=1 时使用
@param EkfBlh 卡尔曼滤波固定解的大地坐标数组，单位度、度、米
@param OutputFile 已打开的绘图数据文件句柄
@param mode 处理模式，1=EKF 2=最小二乘
@return 无
@note 末列状态码 1 表示固定解、2 表示浮点解；两种模式共用一套输出语句，
      仅按模式挑选坐标、Ratio 与固定标志的来源
*/
void WritePlotRecord(SolverData* rawdata, double FixedBlh[3], KalmanState* ekf, double EkfBlh[3], FILE* OutputFile, int mode)
{
    const bool useEkf = (mode == 1);

    const double* blh = useEkf ? EkfBlh : FixedBlh;
    const double ratio = useEkf ? ekf->ratio : rawdata->DDObs.Ratio;
    const bool fixedFlag = useEkf ? ekf->IsFixed : rawdata->DDObs.bFixed;
    const int State = fixedFlag ? 1 : 2;

    fprintf(OutputFile, "%u %.f %0.18f %0.18f %0.18f %0.3f %d\n",
        rawdata->RovEpk_cur.Time.Week, rawdata->RovEpk_cur.Time.SecOfWeek,
        blh[0], blh[1], blh[2],
        ratio, State);
}
