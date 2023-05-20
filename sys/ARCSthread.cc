//! @file ARCSthread.cc
//! @brief ARCSリアルタイムスレッド管理クラス
//!
//! リアルタイムスレッドの生成、開始、停止、破棄などの管理をします。
//!
//! @date 2020/04/15
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <unistd.h>
#include "ARCSthread.hh"
#include "ARCScommon.hh"
#include "ARCSeventlog.hh"
#include "ConstParams.hh"
#include "ARCSassert.hh"
#include "ScreenParams.hh"
#include "GraphPlot.hh"

using namespace ARCS;

//! @brief コンストラクタ
ARCSthread::ARCSthread(ARCSassert& Asrt, ScreenParams& SP, GraphPlot& GP) :
	ARCSast(Asrt),					// ARCSアサートへの参照
	ScrPara(SP),					// 画面パラメータへの参照
	Graph(GP),						// グラフプロットへの参照
	ExpDatMem(),					// 実験データ保存メモリの初期化
	CtrlFuncs(SP, GP, ExpDatMem),	// 制御用周期実行関数群の初期化
	CtrlFuncObj(),					// 制御用周期実行関数の関数オブジェクトの初期化
	RTthreads({nullptr}),			// リアルタイムスレッドへのスマートポインタ配列の初期化
	InfoState(ITS_IDLE),					// 情報取得スレッドの初期状態
	InfoMutex(PTHREAD_MUTEX_INITIALIZER),	// 情報取得スレッド同期用Mutex
	InfoCond(PTHREAD_COND_INITIALIZER),		// 情報取得スレッド同期用条件
	InfoGetThreadID()						// 情報取得スレッドの識別子の初期化
{
	PassedLog();
	pthread_mutex_init(&InfoMutex, nullptr);// 情報取得スレッド同期用Mutex初期化
	pthread_cond_init(&InfoCond, nullptr);	// 情報取得スレッド同期用条件初期化
	CtrlFuncObj = CtrlFuncs.GetCtrlFuncObject();				// 制御用周期実行関数の関数オブジェクトを取得
	for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i){
		RTthreads.at(i) = std::make_unique<SFthread<ConstParams::THREAD_TYPE>>(
			ConstParams::SAMPLING_TIME.at(i), ConstParams::CPUCORE_NUMBER.at(i)
		);	// リアルタイムスレッドの生成
		RTthreads.at(i)->SetRealtimeFunction(CtrlFuncObj[i]);	// 関数オブジェクトをリアルタイムスレッドとして設定
	}
	// 情報取得スレッド生成とCPUコア，ポリシー，優先順位の設定
	pthread_create(&InfoGetThreadID, NULL, (void*(*)(void*))InfoGetThread, this);
	ARCScommon::SetCPUandPolicy(
		InfoGetThreadID,
		ConstParams::ARCS_CPU_INFO,
		ConstParams::ARCS_POL_INFO,
		ConstParams::ARCS_PRIO_INFO
	);
	PassedLog();
}

//! @brief デストラクタ
ARCSthread::~ARCSthread(){
	pthread_mutex_lock(&InfoMutex);		// Mutexロック
	InfoState = ITS_DSTRCT;				// スレッド破棄が指令されたことを知らせる
	pthread_cond_broadcast(&InfoCond);	// 実際の状態が更新されたことをスレッドに知らせる
	pthread_mutex_unlock(&InfoMutex);	// Mutexアンロック
	pthread_join(InfoGetThreadID, nullptr);	// 情報取得スレッド終了待機
	if(ARCSast.IsEmergency() == true){
		// 緊急停止時は，リアルタイムマルチスレッドを強制破壊
		for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->ForceDestruct();
	}
	PassedLog();
}

//! @brief スレッドを開始する関数
void ARCSthread::Start(void){
	// 制御開始！
	PassedLog();
	pthread_mutex_lock(&InfoMutex);		// Mutexロック
	InfoState = ITS_START;				// スレッド開始が指令されたことを知らせる
	pthread_cond_broadcast(&InfoCond);	// 実際の状態が更新されたことをスレッドに知らせる
	pthread_mutex_unlock(&InfoMutex);	// Mutexアンロック
	CtrlFuncs.InitialProcess();	// 初期化モードの実行
	ARCSast.SetRealtimeMode();	// ARCS用assertをリアルタイムモードに変更
	for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->Start();		// リアルタイムマルチスレッドの開始
	for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->WaitStart();	// 開始するまで待機
}

//! @brief スレッドを停止する関数
void ARCSthread::Stop(void){
	// 制御停止！
	EventLog("Waiting for Stop of Realtime Thread...");
	if(ARCSast.IsEmergency() == false){
		// 正常終了時は，
		ARCSast.SetNonRealtimeMode();	// ARCS用assertを非リアルタイムモードに変更
		for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->Stop();	// リアルタイムマルチスレッドの停止
		for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->WaitStop();// リアルタイムマルチスレッドの終了待機
	}else{
		// 緊急停止時は，
		for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->Stop();	// リアルタイムマルチスレッドの停止
		// 終了待機はしない
	}
	EventLog("Waiting for Stop of Realtime Thread...Done");
	CtrlFuncs.ExitProcess();		// 終了処理モードの実行
}

//! @brief スレッドをリセットする関数
void ARCSthread::Reset(void){
	for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i) RTthreads.at(i)->Reset();	// リアルタイムマルチスレッドのリセット
	ExpDatMem.Reset();	// 実験データ保存メモリもリセット
}

//! @brief 測定データを保存する関数
void ARCSthread::SaveDataFiles(void){
	EventLog("Writing PNG/CSV Data Files...");
	Graph.SaveScreenImage();	// スクリーンショットをPNGに保存
	ExpDatMem.WriteCsvFile();	// データメモリの中身をCSVに保存
	EventLog("Writing PNG/CSV Data Files...Done");
}

//! @brief 情報取得スレッド
//! @param[in]	p	クラスメンバアクセス用ポインタ
void ARCSthread::InfoGetThread(ARCSthread* const p){
	double Time = 0;													// [s] 時刻 (一番速いスレッド RTthreads[0] の時刻)
	std::array<double, ConstParams::THREAD_MAX> PeriodicTime = {0};		// [s] 計測された制御周期
	std::array<double, ConstParams::THREAD_MAX> ComputationTime = {0};	// [s] 計測された消費時間
	std::array<double, ConstParams::THREAD_MAX> MaxTime = {0};			// [s] 計測された制御周期の最大値
	std::array<double, ConstParams::THREAD_MAX> MinTime = {0};			// [s] 計測された制御周期の最小値
	
	// 動作状態が「開始」か「破棄」に設定されるまで待機
	EventLog("Waiting for ITS_START,ITS_DSTRCT...");
	pthread_mutex_lock(&(p->InfoMutex));	// Mutexロック
	while(p->InfoState != ITS_START && p->InfoState != ITS_DSTRCT){
		pthread_cond_wait(&(p->InfoCond), &(p->InfoMutex));	// 状態が更新されるまで待機
	}
	pthread_mutex_unlock(&(p->InfoMutex));	// Mutexアンロック
	EventLog("Waiting for ITS_START,ITS_DSTRCT...Done");
	
	if(p->InfoState == ITS_DSTRCT) return;	// 破棄指令ならスレッド終了
	
	// スレッドが破棄されるまで無限ループ
	while(p->InfoState != ITS_DSTRCT){
		// リアルタイムスレッドから時刻情報を取得
		Time = p->RTthreads.at(0)->GetTime();	// 時刻の取得
		
		for(size_t i = 0; i < ConstParams::THREAD_NUM; ++i){
			PeriodicTime[i]    = p->RTthreads[i]->GetSmplTime();// 制御周期の取得
			ComputationTime[i] = p->RTthreads[i]->GetCompTime();// 消費時間の取得
			MaxTime[i]         = p->RTthreads[i]->GetMaxTime();	// 制御周期の最大値の取得
			MinTime[i]         = p->RTthreads[i]->GetMinTime();	// 制御周期の最小値の取得
		}
		
		// ARCS画面パラメータに格納
		p->ScrPara.SetTime(Time);
		p->ScrPara.SetTimeVars(PeriodicTime, ComputationTime, MaxTime, MinTime);
		
		// 制御用変数値を更新
		p->CtrlFuncs.UpdateControlValue();	// 制御用周期実行関数群クラス内部で，ARCS画面パラメータに書き込む/から読み込む
		
		usleep(ConstParams::ARCS_TIME_INFO);// 指定時間だけ待機
	}
	
	EventLog("InfoGetThread Destructed.");
}

