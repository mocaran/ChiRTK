#include "gnss_core_types.h"


/*!
@brief 采用Hopfield模型计算对流层斜路径延迟改正
@param H 测站高程（米），相对大地水准面
@param Elev 卫星高度角（度）
@return 对流层延迟改正量（米），测站高于对流层顶时返回0
@note 干分量与湿分量分别按不同映射函数计算后求和
*/
double HopfieldTroposphereDelay(const double H, const double Elev)
{
	// 容错处理：如果测站高于对流层顶，返回0
	if (abs(H) > Htrop) {
		return 0.0;
	}

	// 计算准备参数
	double RH = RelHum * exp(-0.0006396 * (H - SeaLevel));
	double p = AtmosPre * pow(1.0 - 0.0000226 * (H - SeaLevel), 5.225);

	double T = TempT0 - 0.0065 * (H - SeaLevel);
	double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
	double h_w = 11000.0;
	double h_d = 40136.0 + 148.72 * (TempT0 - 273.16);
	double K_w = 155.2e-7 * 4810.0 * e * (h_w - H) / pow(T, 2);
	double K_d = 155.2e-7 * p * (h_d - H) / T;

	// 计算delta
	double delta_d, delta_w;
	double Sin_1 = sin(sqrt(pow(Elev, 2) + 6.25) * Rad);
	double Sin_2 = sin(sqrt(pow(Elev, 2) + 2.25) * Rad);
	delta_d = K_d / Sin_1;
	delta_w = K_w / Sin_2;

	// 返回对流层改正值
	return delta_d + delta_w;
}



/*!
@brief 用GF与MW组合的历元间差值探测并标记粗差
@param Obs 当前历元观测数据，读写SatObs.Valid与ComObs
@return 无
@note GF与MW阈值分别为5厘米和3米，超限卫星标记为无效
*/
void RejectObservationOutliers(EpochObservation* Obs)
{
	// 当前历元GF和MW临时结算结果
	CombinationState CurComObs[MAXCHANNUM];
	memset(CurComObs, 0, sizeof(CurComObs));

	for (int i = 0; i < Obs->SatNum; i++)
	{
		Obs->SatObs[i].Valid = false;

		// 1、检查该卫星的双频伪距和相位数据是否有效和完整，若不全或为0，将Valid标记为false，continue
		if (fabs(Obs->SatObs[i].P[0]) < 1e-5 || fabs(Obs->SatObs[i].P[1]) < 1e-5 || fabs(Obs->SatObs[i].L[0]) < 1e-5 || fabs(Obs->SatObs[i].L[1]) < 1e-5)
		{
			continue;
		}

		// 2、计算当前历元该卫星的GF和MW组合值
		double Lgf, Lmw;

		if (Obs->SatObs[i].System == GPS)	// 如果是GPS系统
		{
			Lgf = Obs->SatObs[i].L[0] - Obs->SatObs[i].L[1];
			Lmw = (FG1_GPS * Obs->SatObs[i].L[0] - FG2_GPS * Obs->SatObs[i].L[1]) / (FG1_GPS - FG2_GPS)
				- (FG1_GPS * Obs->SatObs[i].P[0] + FG2_GPS * Obs->SatObs[i].P[1]) / (FG1_GPS + FG2_GPS);

		}
		else if (Obs->SatObs[i].System == BDS)	// 如果是BDS系统
		{
			Lgf = Obs->SatObs[i].L[0] - Obs->SatObs[i].L[1];
			Lmw = (FG1_BDS * Obs->SatObs[i].L[0] - FG3_BDS * Obs->SatObs[i].L[1]) / (FG1_BDS - FG3_BDS)
				- (FG1_BDS * Obs->SatObs[i].P[0] + FG3_BDS * Obs->SatObs[i].P[1]) / (FG1_BDS + FG3_BDS);
		}
		else
		{
			Obs->SatObs[i].Valid = false;
			continue;
		}

		// 存储当前MWGF计算值
		CurComObs[i].Prn = Obs->SatObs[i].Prn;
		CurComObs[i].Sys = Obs->SatObs[i].System;
		CurComObs[i].GF = Lgf;
		CurComObs[i].MW = Lmw;
		CurComObs[i].n = 1;


		// 3、对于可用的观测数据，计算伪距的IF组合观测值，用于SPP
		if (CurComObs[i].Sys == GPS)
		{
			CurComObs[i].PIF = (FG1_GPS * FG1_GPS * Obs->SatObs[i].P[0] - FG2_GPS * FG2_GPS * Obs->SatObs[i].P[1]) / (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS);
		}
		else if (CurComObs[i].Sys == BDS)
		{
			CurComObs[i].PIF = (FG1_BDS * FG1_BDS * Obs->SatObs[i].P[0] - FG3_BDS * FG3_BDS * Obs->SatObs[i].P[1]) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);
		}

		CurComObs[i].n = 1;	// 作为新历史的初始计数
		for (int j = 0; j < MAXCHANNUM; j++)
		{
			// 4、从上个历元的MWGF数据中查找该卫星的GF和MW组合值
			if (Obs->ComObs[j].Prn == CurComObs[i].Prn && Obs->ComObs[j].Sys == CurComObs[i].Sys)
			{
				// 5、6、计算当前历元该卫星GF与上一历元对应GF的差值dGF，计算当前历元该卫星MW与上一历元对应MW平滑值的差值dMW
				double dGF = fabs(CurComObs[i].GF - Obs->ComObs[j].GF);
				double dMW = fabs(CurComObs[i].MW - Obs->ComObs[j].MW);

				// 7、检查dGF和dMW是否超限，限差阈值建议为5cm和3m。若超限，标记为粗差，将Valid标记为false ，若不超限，标记为可用将Valid标记为true，并计算该卫星的MW平滑值
				if (dGF < GFThres && dMW < MWThres)
				{
					Obs->SatObs[i].Valid = true;
					CurComObs[i].MW = (Obs->ComObs[j].MW * Obs->ComObs[j].n + CurComObs[i].MW) / (Obs->ComObs[j].n + 1);
					CurComObs[i].n = Obs->ComObs[j].n + 1;
				}
				// 本身初始化的时候就是false，不需要处理
				break;
			}
		}
	}
	// 8、将缓冲区组合观测值的拷贝到obs里
	memcpy(Obs->ComObs, CurComObs, sizeof(CurComObs));
}
