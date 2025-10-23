//! @file ARCSscrparams.cc
//! @brief ARCS画面パラメータ格納クラス
//!        ARCS用画面に表示する各種パラメータを格納します。
//! @date 2024/06/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCS_MR

#include <tuple>
#include "ARCSscrparams.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ
ARCSscrparams::ARCSscrparams(void)
	: TimeMutex(PTHREAD_MUTEX_INITIALIZER),
	  ActMutex(PTHREAD_MUTEX_INITIALIZER),
	  IndicMutex(PTHREAD_MUTEX_INITIALIZER),
	  OnsetMutex(PTHREAD_MUTEX_INITIALIZER),
	  Time(0),
	  PeriodicTime({0}),
	  ComputationTime({0}),
	  MaxTime({0}),
	  MinTime({0}),
	  NetworkLink(false),
	  Initializing(false),
	  CurrentRef({0}),
	  PositionRes({0}),
	  VarIndicator({0}),
	  VarIndicatorBuf({0}),
	  VarIndicCount(0),
	  OnlineSetVar({0}),
	  OnlineSetVarIni({0}),
	  SetVarCount(0),
	  InitSetVarCount(0)
{
	PassedLog();
	// Mutexの初期化
	pthread_mutex_init(&TimeMutex, nullptr);
	pthread_mutex_init(&ActMutex, nullptr);
	pthread_mutex_init(&IndicMutex, nullptr);
	pthread_mutex_init(&OnsetMutex, nullptr);
}

//! @brief デストラクタ
ARCSscrparams::~ARCSscrparams(){
	PassedLog();
}

//! @brief 時刻を取得する関数
//! @return	時刻
double ARCSscrparams::GetTime(void){
	double ret;
	pthread_mutex_lock(&TimeMutex);
	ret = Time;
	pthread_mutex_unlock(&TimeMutex);
	return ret;
}

//! @brief 時刻を設定する関数
//! @param[in]	t	時刻
void ARCSscrparams::SetTime(const double t){
	pthread_mutex_lock(&TimeMutex);
	Time = t;
	pthread_mutex_unlock(&TimeMutex);
}

//! @brief 実際の制御周期，消費時間，制御周期の最大値，最小値を返す関数
//! @param[in]	ThreadNum	リアルタイムスレッド番号
//! @return 制御周期，消費時間，制御周期の最大値，最小値
std::tuple<double, double, double, double> ARCSscrparams::GetTimeVars(const unsigned int ThreadNum){
	double PT, CT, Max, Min;
	pthread_mutex_lock(&TimeMutex);
	PT  = PeriodicTime.at(ThreadNum);
	CT  = ComputationTime.at(ThreadNum);
	Max = MaxTime.at(ThreadNum);
	Min = MinTime.at(ThreadNum);
	pthread_mutex_unlock(&TimeMutex);
	return std::forward_as_tuple(PT, CT, Max, Min);
}

//! @brief 実際の制御周期，消費時間，制御周期の最大値，最小値の配列を返す関数
//! @param[out]	PT	制御周期の配列
//! @param[out]	CT	消費時間の配列
//! @param[out]	Max	制御周期の最大値の配列
//! @param[out]	Min	制御周期の最小値の配列
void ARCSscrparams::GetTimeVars(
	std::array<double, ARCSparams::THREAD_MAX>& PT,
	std::array<double, ARCSparams::THREAD_MAX>& CT,
	std::array<double, ARCSparams::THREAD_MAX>& Max,
	std::array<double, ARCSparams::THREAD_MAX>& Min
){
	pthread_mutex_lock(&TimeMutex);
	PT  = PeriodicTime;
	CT  = ComputationTime;
	Max = MaxTime;
	Min = MinTime;
	pthread_mutex_unlock(&TimeMutex);
}

//! @brief 実際の制御周期，消費時間，制御周期の最大値，最小値の配列を設定する関数
//! @param[in]	PT	制御周期の配列
//! @param[in]	CT	消費時間の配列
//! @param[in]	Max	制御周期の最大値の配列
//! @param[in]	Min	制御周期の最小値の配列
void ARCSscrparams::SetTimeVars(
	const std::array<double, ARCSparams::THREAD_MAX>& PT,
	const std::array<double, ARCSparams::THREAD_MAX>& CT,
	const std::array<double, ARCSparams::THREAD_MAX>& Max,
	const std::array<double, ARCSparams::THREAD_MAX>& Min
){
	pthread_mutex_lock(&TimeMutex);
	PeriodicTime = PT;
	ComputationTime = CT;
	MaxTime = Max;
	MinTime = Min;
	pthread_mutex_unlock(&TimeMutex);
}

//! @brief ネットワークリンクフラグを取得する関数
//! @return ネットワークリンクフラグ
bool ARCSscrparams::GetNetworkLink(void){
	return NetworkLink;
}

//! @brief ネットワークリンクフラグを設定する関数
//! @param[in] LinkFlag	ネットワークリンクフラグ
void ARCSscrparams::SetNetworkLink(const bool LinkFlag){
	NetworkLink = LinkFlag;
}

//! @brief ロボット初期化フラグを取得する関数
//! @return ロボット初期化フラグ
bool ARCSscrparams::GetInitializing(void){
	return Initializing;
}

//! @brief ロボット初期化フラグを設定する関数
//! @param[in] InitFlag	ロボット初期化フラグ
void ARCSscrparams::SetInitializing(const bool InitFlag){
	Initializing = InitFlag;
}

//! @brief 電流と位置を取得する関数
//! @param[in]	ActNum	アクチュエータ番号
//! @return	電流指令，位置応答
std::tuple<double, double> ARCSscrparams::GetCurrentAndPosition(const unsigned int ActNum){
	double Current, Position;
	pthread_mutex_lock(&ActMutex);
	Current  = CurrentRef.at(ActNum);
	Position = PositionRes.at(ActNum);
	pthread_mutex_unlock(&ActMutex);
	return std::forward_as_tuple(Current, Position);
}

//! @brief 電流と位置の配列を取得する関数
//! @param[in]	Current	電流指令の配列
//! @param[in]	Position	位置応答の配列
void ARCSscrparams::GetCurrentAndPosition(
	std::array<double, EquipParams::ACTUATOR_NUM>& Current,
	std::array<double, EquipParams::ACTUATOR_NUM>& Position
){
	pthread_mutex_lock(&ActMutex);
	Current  = CurrentRef;
	Position = PositionRes;
	pthread_mutex_unlock(&ActMutex);
}

//! @brief 電流と位置の配列を設定する関数
//! @param[in]	Current		電流指令ベクトル
//! @param[in]	Position	位置ベクトル
void ARCSscrparams::SetCurrentAndPosition(
	const ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Current,
	const ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Position
){
	pthread_mutex_lock(&ActMutex);
	Current.StoreArray(CurrentRef);
	Position.StoreArray(PositionRes);
	pthread_mutex_unlock(&ActMutex);
}

//! @brief 任意変数インジケータの配列を返す関数
//! @param[out]	Vars	任意変数値の配列
void ARCSscrparams::GetVarIndicator(std::array<double, ARCSparams::INDICVARS_MAX>& Vars){
	pthread_mutex_lock(&IndicMutex);
	Vars = VarIndicator;
	pthread_mutex_unlock(&IndicMutex);
}

//! @brief 任意変数インジケータの配列を設定する関数
//! @param[in]	Vars	任意変数値の配列
void ARCSscrparams::SetVarIndicator(const std::array<double, ARCSparams::INDICVARS_MAX>& Vars){
	pthread_mutex_lock(&IndicMutex);
	VarIndicator = Vars;
	pthread_mutex_unlock(&IndicMutex);
}

//! @brief オンライン設定変数の配列を返す関数
//! @param[out]	Vars	オンライン設定変数値の配列
void ARCSscrparams::GetOnlineSetVars(std::array<double, ARCSparams::ONLINEVARS_MAX>& Vars){
	pthread_mutex_lock(&OnsetMutex);
	Vars = OnlineSetVar;
	pthread_mutex_unlock(&OnsetMutex);
}

//! @brief オンライン設定変数の配列を設定する関数
//! @param[in]	Vars	オンライン設定変数値の配列
void ARCSscrparams::SetOnlineSetVars(const std::array<double, ARCSparams::ONLINEVARS_MAX>& Vars){
	pthread_mutex_lock(&OnsetMutex);
	OnlineSetVar = Vars;
	pthread_mutex_unlock(&OnsetMutex);
}

//! @brief オンライン設定変数に値を設定する関数
//! @param[in]	VarNum	変数番号
//! @param[in]	VarVal	オンライン設定変数値
void ARCSscrparams::SetOnlineSetVar(const unsigned int VarNum, const double VarVal){
	pthread_mutex_lock(&OnsetMutex);
	OnlineSetVar[VarNum] = VarVal;
	pthread_mutex_unlock(&OnsetMutex);
}

#else // ARCS_MR

#include "ARCSscrparams.hh"
#include "WebXR.hh"

using namespace ARCS;

void ARCSscrparams::SetNetworkLink(const bool LinkFlag) {}

void ARCSscrparams::SetInitializing(const bool InitFlag) {}

void ARCSscrparams::SetCurrentAndPosition(
		const ArcsMat<ConstParams::ACTUATOR_NUM, 1>& Current,
		const ArcsMat<ConstParams::ACTUATOR_NUM, 1>& Position
	) {
  for (int i = 1; i <= ConstParams::ACTUATOR_NUM; ++i) {
    int status = ACTUATOR_STATUS_NORMAL;

    setActuatorStatus(i, status, Current(i, 1), Position(i, 1));
  }
}

void ARCSscrparams::SetVarIndicator(
    const std::array<double, ConstParams::INDICVARS_MAX> &Vars) {
  for (int i = 0; i < ConstParams::INDICVARS_NUM; ++i) {
    setReadVariable(i, Vars[i]);
  }
}

void ARCSscrparams::GetOnlineSetVars(
    std::array<double, ConstParams::ONLINEVARS_MAX> &Vars) {
  for (int i = 0; i < ConstParams::INDICVARS_NUM; ++i) {
    if (canGetWriteVariable(i)) {
      Vars[i] = getWriteVariable(i);
    }
  }
}

void ARCSscrparams::SetOnlineSetVars(
    const std::array<double, ConstParams::ONLINEVARS_MAX> &Vars) {
  for (int i = 0; i < ConstParams::ONLINEVARS_NUM; ++i) {
    setWriteVariable(i, Vars[i]);
  }
}

void ARCSscrparams::SetOnlineSetVar(const unsigned int VarNum,
                                   const double VarVal) {
  setWriteVariable(VarNum, VarVal);
}

void ARCSscrparams::UpdateOnlineSetVar() { GetOnlineSetVars(OnlineSetVar); }

#endif // ARCS_MR