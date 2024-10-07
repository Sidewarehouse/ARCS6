//! @file SFthread.hh
//! @brief SCHED_FIFOリアルタイムスレッドクラス(sleep使用不使用テンプレート可変版, 関数オブジェクト版)
//!
//! pthreadのSCHED_FIFOで実時間スレッドを生成＆管理＆破棄する。実際に計測された制御周期や計算消費時間も提供する。
//!
//! @date 2024/10/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef SFTHREADING
#define SFTHREADING

#include <sys/mman.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <functional>
#include <cfenv>
#include <cstdint>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "CPUSettings.hh"
#include "LinuxCommander.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

namespace ARCS {	// ARCS名前空間
//! @brief リアルタイムアルゴリズムのタイプの定義
//! 解説：
//! 時間ゼロのスリープを入れるとLinuxが安定動作できるが，リアルタイム性は犠牲になる。
//! 古い低速の処理系では INSERT_ZEROSLEEP を推奨する。
//! 新しい高速の処理系では WITHOUT_ZEROSLEEP に設定するとリアルタイム性が改善する。
enum class SFalgorithm {
	INSERT_ZEROSLEEP,	//!< リアルタイムループに時間ゼロのスリープを入れる
	WITHOUT_ZEROSLEEP	//!< スリープは一切入れない
};

//! @brief カーネルパラメータのタイプの定義
enum class SFkernelparam : uint8_t {
	NO_SETTINGS     = 0b00000000,	//!< 何も設定しない場合
	CFS_DISABLED	= 0b00000001,	//!< CFS(Completely Fair Scheduler)を無効にする
	PREEMPT_DYNFULL = 0b00000010	//!< PREEMPT_DYNAMICの場合にFULLモードにする
};
// ↑を併用する場合は、
// static_cast<SFkernelparam>( static_cast<uint8_t>(SFkernelparam::CFS_DISABLED) | static_cast<uint8_t>(SFkernelparam::PREEMPT_DYNFULL) )
// のようにキャストを駆使して論理和を取ること。

//! @brief 実時間スレッド生成・破棄クラス
//! @tparam	SFA	リアルタイムアルゴリズムのタイプ
template <SFalgorithm SFA = SFalgorithm::INSERT_ZEROSLEEP, SFkernelparam SFK = SFkernelparam::CFS_DISABLED>
class SFthread {
	public:
		//! @brief 動作状態の定義
		enum ThreadState {
			SFID_ERROR,	//!< エラー検出
			SFID_START,	//!< 開始
			SFID_RUN,	//!< 動作中
			SFID_STOP,	//!< 停止
			SFID_EXCMPL,//!< 終了動作完了
			SFID_DSTRCT	//!< スレッド破棄
		};
		
		//! @brief コンストラクタ
		//! @param[in] PeriodTime	制御周期
		//! @param[in] FuncObject	制御用実行関数の関数オブジェクト
		//! @param[in] CPUno		使用するCPUコアの番号
		SFthread(const unsigned long PeriodTime, const std::function<bool(double,double,double)>& FuncObject, const int CPUno)
			: SyncMutex(PTHREAD_MUTEX_INITIALIZER),	// 同期用Mutex
			  SyncCond(PTHREAD_COND_INITIALIZER),	// 同期用条件
			  StateFlag(SFID_STOP),		// 動作状態フラグを「停止状態」に設定
			  Ts(PeriodTime),			// [ns] 制御周期の格納
			  FuncObj(FuncObject),		// 制御用実行関数への関数オブジェクトを格納
			  Time(),					// 時刻の初期化
			  ActPeriodicTime(),		// 周期時間の初期化
			  ComputationTime(),		// 消費時間の初期化
			  ThreadID(0),				// スレッド識別子の初期化
			  ThreadParam(),			// スレッドパラメータ
			  MaxMemo(0),				// サンプリング時間最大値計算用
			  MinMemo(PeriodTime*1e-9)	// サンプリング時間最小値計算用
		{
			// 実時間スレッドの生成と優先度の設定
			PassedLog();
			if constexpr(SFA == SFalgorithm::INSERT_ZEROSLEEP){
				EventLog("INSERT_ZEROSLEEP MODE.");
			}else{
				EventLog("WITHOUT_ZEROSLEEP MODE.");
			}
			SetKernelParameters();		// カーネルパラメータをリアルタイム用に設定
			pthread_mutex_init(&SyncMutex, nullptr);	// 同期用Mutexの初期化
			pthread_cond_init(&SyncCond, nullptr);		// 同期用条件の初期化
			pthread_create(&ThreadID, NULL, (void*(*)(void*))RealTimeThread, this);	// スレッド生成
			CPUSettings::SetCPUandPolicy(ThreadID, CPUno, SCHED_FIFO);				// CPUコアの割り当てとスケジューリングポリシーの設定
			PassedLog();
		}
		
		//! @brief コンストラクタ(制御用実行関数の関数オブジェクトを後で設定する場合)
		//! @param[in] PeriodTime	制御周期
		//! @param[in] CPUno		使用するCPUコアの番号
		SFthread(const unsigned long PeriodTime, const int CPUno)
			: SyncMutex(PTHREAD_MUTEX_INITIALIZER),	// 同期用Mutex
			  SyncCond(PTHREAD_COND_INITIALIZER),	// 同期用条件
			  StateFlag(SFID_STOP),		// 動作状態フラグを「停止状態」に設定
			  Ts(PeriodTime),			// [ns] 制御周期の格納
			  FuncObj(),				// 制御用実行関数への関数オブジェクトを格納
			  Time(),					// 時刻の初期化
			  ActPeriodicTime(),		// 周期時間の初期化
			  ComputationTime(),		// 消費時間の初期化
			  ThreadID(0),				// スレッド識別子の初期化
			  ThreadParam(),			// スレッドパラメータ
			  MaxMemo(0),				// サンプリング時間最大値計算用
			  MinMemo(PeriodTime*1e-9)	// サンプリング時間最小値計算用
		{
			// 実時間スレッドの生成と優先度の設定
			PassedLog();
			if constexpr(SFA == SFalgorithm::INSERT_ZEROSLEEP){
				EventLog("INSERT_ZEROSLEEP MODE.");
			}else{
				EventLog("WITHOUT_ZEROSLEEP MODE.");
			}
			SetKernelParameters();		// カーネルパラメータをリアルタイム用に設定
			pthread_mutex_init(&SyncMutex, nullptr);	// 同期用Mutexの初期化
			pthread_cond_init(&SyncCond, nullptr);		// 同期用条件の初期化
			pthread_create(&ThreadID, NULL, (void*(*)(void*))RealTimeThread, this);	// スレッド生成
			CPUSettings::SetCPUandPolicy(ThreadID, CPUno, SCHED_FIFO);				// CPUコアの割り当てとスケジューリングポリシーの設定
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r 右辺値
		SFthread(SFthread&& r)
			: SyncMutex(r.SyncMutex),			// 同期用Mutex
			SyncCond(r.SyncCond),				// 同期用条件
			StateFlag(r.StateFlag),				// 動作状態フラグ
			Ts(r.Ts),							// [ns] 制御周期
			FuncObj(r.FuncObj),					// 制御用実行関数への関数オブジェクト
			Time(r.Time),						// 時刻
			ActPeriodicTime(r.ActPeriodicTime),	// 周期時間
			ComputationTime(r.ComputationTime),	// 消費時間
			ThreadID(r.ThreadID),				// スレッド識別子
			ThreadParam(r.ThreadParam),			// スレッドパラメータ
			MaxMemo(r.MaxMemo),					// サンプリング時間最大値計算用
			MinMemo(r.MinMemo)					// サンプリング時間最小値計算用
		{
			
		}
		
		//! @brief デストラクタ
		~SFthread(){
			PassedLog();
			pthread_mutex_lock(&SyncMutex);		// Mutexロック
			StateFlag = SFID_DSTRCT;			// スレッド破棄が指令されたことを知らせる
			pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことをリアルタイムスレッドに知らせる
			pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
			pthread_join(ThreadID, nullptr);	// 実時間スレッド終了待機
			ResetKernelParameters();			// カーネルパラメータを元に戻す
			PassedLog();
		}
		
		//! @brief 実時間スレッドから呼び出す関数を設定する関数
		//! @param[in]	FuncObject	関数オブジェクト
		void SetRealtimeFunction(std::function<bool(double,double,double)> FuncObject){
			FuncObj = FuncObject;	// 関数オブジェクトをセット
		}
		
		//! @brief スレッド実行を開始する関数
		void Start(void){
			pthread_mutex_lock(&SyncMutex);		// Mutexロック
			StateFlag = SFID_START;				// 開始が指令されたことを知らせる
			pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことをリアルタイムスレッドに知らせる
			pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
		}
		
		//! @brief スレッド実行が開始されるまで待機する関数
		void WaitStart(void){
			// 状態が「動作中」に設定されるまで待機
			EventLog("Waiting for SFID_RUN...");
			pthread_mutex_lock(&SyncMutex);		// Mutexロック
			while(StateFlag != SFID_RUN){
				pthread_cond_wait(&SyncCond, &SyncMutex);	// 状態が更新されるまで待機
			}
			pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
			EventLog("Waiting for SFID_RUN...Done");
		}
		
		//! @brief スレッド実行を停止する関数
		void Stop(void){
			pthread_mutex_lock(&SyncMutex);		// Mutexロック
			StateFlag = SFID_STOP;				// 停止が指令されたことを知らせる
			pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことをリアルタイムスレッドに知らせる
			pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
		}
		
		//! @brief スレッド実行が停止されるまで待機する関数
		void WaitStop(void){
			// 状態が「終了動作完了」に設定されるまで待機
			EventLog("Waiting for SFID_EXCMPL...");
			pthread_mutex_lock(&SyncMutex);		// Mutexロック
			while(StateFlag != SFID_EXCMPL){
				pthread_cond_wait(&SyncCond, &SyncMutex);	// 状態が更新されるまで待機
			}
			pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
			EventLog("Waiting for SFID_EXCMPL...Done");
		}
		
		//! @brief スレッドをリセットする関数
		void Reset(void){
			timespec_clear(Time);			// 時刻をクリア
			timespec_clear(ActPeriodicTime);// 実際の周期時間をクリア
			timespec_clear(ComputationTime);// 消費時間をクリア
			MaxMemo = 0;		// 計測周期最大値をクリア
			MinMemo = Ts*1e-9;	// 計測周期最小値をクリア
		}
		
		//! @brief スレッドを強制破壊する関数
		void ForceDestruct(void){
			pthread_cancel(ThreadID);		// スレッド破棄
			pthread_join(ThreadID, nullptr);// 実時間スレッド終了待機
		}
		
		//! @brief 時刻を取得する関数
		//! @return 時刻 [s]
		double GetTime(void) const {
			return timespec_to_sec(Time);			// [s] 時刻に変換して返す
		}
		
		//! @brief 計測された実際のサンプリング時間を取得する関数
		//! @return 計測周期 [s]
		double GetSmplTime(void) const {
			return timespec_to_sec(ActPeriodicTime);// [s] 時刻に変換して返す
		}
		
		//! @brief 計測された消費時間を取得する関数
		//! @return 計測消費時間 [s]
		double GetCompTime(void) const {
			return timespec_to_sec(ComputationTime);// [s] 時刻に変換して返す
		}
		
		//! @brief 計測された実際のサンプリング時間の最大値を取得する関数
		//! @return 計測最大サンプリング時間 [s]
		double GetMaxTime(void){
			double TsZ0 = GetSmplTime();	// [s] 今のサンプリング時間を取得
			if(MaxMemo < TsZ0){
				MaxMemo = TsZ0;
			}
			return MaxMemo;
		}
		
		//! @brief 計測された実際のサンプリング時間の最小値を取得する関数
		//! @return 計測最小サンプリング時間 [s]
		double GetMinTime(void){
			double TsZ0 = GetSmplTime();	// [s] 今のサンプリング時間を取得
			if(TsZ0 < MinMemo && 1e-6 < TsZ0){
				MinMemo = TsZ0;
			}
			return MinMemo;
		}
		
	private:
		SFthread(const SFthread&) = delete;					//!< コピーコンストラクタ使用禁止
		const SFthread& operator=(const SFthread&) = delete;//!< 代入演算子使用禁止
		
		static const long ONE_SEC_IN_NANO = 1000000000;		//!< [ns] 1秒をナノ秒で表すと
		static const long SKEW_IN_NANO = 40;				//!< [ns] 時刻ズレ調整用パラメータ
		pthread_mutex_t SyncMutex;							//!< 同期用Mutex
		pthread_cond_t	SyncCond;							//!< 同期用条件
		enum ThreadState StateFlag;							//!< 動作状態フラグ
		const unsigned long Ts;								//!< 制御周期
		std::function<bool(double,double,double)> FuncObj;	//!< 関数オブジェクト 引数(時刻, 計測周期, 消費時間)
		timespec Time;										//!< 計測された実際の時刻
		timespec ActPeriodicTime;							//!< 計測された実際の周期時間
		timespec ComputationTime;							//!< 計算によって消費された時間 (つまり ComputationTime < ActPeriodicTime でなければならない)
		pthread_t ThreadID;									//!< スレッド識別子
		struct sched_param ThreadParam;						//!< スレッドパラメータ
		double MaxMemo;										//!< [s] サンプリング時間最大値計算用
		double MinMemo;										//!< [s] サンプリング時間最小値計算用
		
		//! @brief リアルタイムループ
		//! 実際の制御用実行関数はこの関数から呼ばれている
		void RealTimeLoop(void){
			timespec InitTime = {0};		// ループに入る初期時刻
			timespec NextTime = {0};		// 次のループ開始時刻
			timespec TimeInWait = {0};		// 待機ループ内で取得する現在時刻
			timespec PeriodTime = nsec_to_timespec(Ts - SKEW_IN_NANO);	// 所望の制御周期(時刻ズレ調整済み)
			timespec StartTime = {0};		// 開始時刻格納用
			timespec StartTimePrev = {0};	// 前回の開始時間格納用
			timespec EndTime = {0};			// 終了時刻格納用
			timespec PreventStuck = {0};	// 「BUG: soft lockup - CPU#0 Stuck for 67s!」を回避するためのスリープ用
			bool ClockOverride = false;		// 時刻待機のクロックオーバーライドフラグ
			
			EventLog("Starting Realtime Loop.");
			
			clock_gettime(CLOCK_MONOTONIC, &InitTime);			// 初期開始時刻の取得
			StartTimePrev = timespec_sub(InitTime, PeriodTime);	// 実際の制御周期計算用の初期値設定
			
			// 実時間ループ
			while(StateFlag != SFID_STOP){	// 動作状態フラグが「停止」に設定されるまでループ
				// ここからリアルタイム空間
				clock_gettime(CLOCK_MONOTONIC, &StartTime);							// 開始時刻の取得
				Time = timespec_sub(StartTime, InitTime);							// 実際の時刻を計算
				ActPeriodicTime = timespec_sub(StartTime, StartTimePrev);			// 実際の周期時間を計算(timespec構造体は単純に減算できないことに注意)

				std::feclearexcept(FE_ALL_EXCEPT);									// 浮動小数点例外フラグをクリア
				ClockOverride = !FuncObj(GetTime(), GetSmplTime(), GetCompTime());	// 制御用関数の実行(関数オブジェクトにより、ここで実際の制御関数が呼ばれる)
				arcs_assert(std::fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW) == false);	// 浮動小数点例外チェック(ゼロ割、NaN、桁溢れ検出)
				StartTimePrev = StartTime;											// 次回用に今回の開始時刻を格納
				NextTime = timespec_add(StartTime, PeriodTime);						// 開始時刻に制御周期を加算して次の時刻を計算
				
				if constexpr(SFA == SFalgorithm::INSERT_ZEROSLEEP){
					clock_nanosleep(CLOCK_MONOTONIC, 0, &PreventStuck, nullptr);	// 「BUG: soft lockup - CPU#0 Stuck for 67s!」を回避するためのスリープ
				}
				clock_gettime(CLOCK_MONOTONIC, &EndTime);							// 終了時刻の取得
				ComputationTime = timespec_sub(EndTime, StartTime);					// 消費時間を計算(timespec構造体は単純に減算できないことに注意)
				
				// 次の時刻になるまで待機
				while(StateFlag != SFID_STOP){
					clock_gettime(CLOCK_MONOTONIC, &TimeInWait);					// 現在時刻の取得
					if(timespec_lessthaneq(NextTime, TimeInWait) == true || ClockOverride == true){
						// 現在時刻と予め計算した次の時刻とを比較して，超えたら待機終了
						// もしくは，クロックオーバーライドフラグが立っていたら待機終了
						break;
					}
				}
				// リアルタイム空間ここまで
			}	
			
			EventLog("Ending Realtime Loop.");
		}
		
		//! @brief リアルタイムスレッド
		//! @param[in]	p	クラスメンバアクセス用ポインタ
		static void RealTimeThread(SFthread *p){
			// デストラクタが呼ばれるまで繰り返し続ける
			while(1){
				// 動作状態が「開始」か「破棄」に設定されるまで待機
				EventLog("Waiting for SFID_START,SFID_DSTRCT...");
				pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
				while(p->StateFlag != SFID_START && p->StateFlag != SFID_DSTRCT){
					pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 状態が更新されるまで待機
				}
				pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
				EventLog("Waiting for SFID_START,SFID_DSTRCT...Done");
				
				if(p->StateFlag == SFID_DSTRCT) break;	// 破棄指令ならスレッド終了
				
				pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
				p->StateFlag = SFID_RUN;				// 動作状態フラグを「動作中」に設定
				pthread_cond_broadcast(&(p->SyncCond));	// 実際の状態が更新されたことを上位系に知らせる
				pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
				
				mlockall(MCL_CURRENT);					// メモリロック
				p->RealTimeLoop();						// 実時間ループの実行
				munlockall();							// メモリロック解除
				
				pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
				p->StateFlag = SFID_EXCMPL;				// 動作状態フラグを「終了動作完了」に設定
				pthread_cond_broadcast(&(p->SyncCond));	// 実際の状態が更新されたことを上位系に知らせる
				pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
			}
			EventLog("RealTimeThread Destructed.");
		}
		
		//! @brief timespec構造体用の時刻を加算する関数 a + b
		//! @param[in]	a	timespec時刻
		//! @param[in]	b	timespec時刻
		//! @return 加算結果
		static timespec timespec_add(const timespec& a, const timespec& b){
			timespec ret;
			ret.tv_sec = a.tv_sec + b.tv_sec;		// aとbの秒カウンタを加算
			ret.tv_nsec = a.tv_nsec + b.tv_nsec;	// aとbのナノ秒カウンタを加算
			if(ONE_SEC_IN_NANO <= ret.tv_nsec){
				// ナノ秒カウンタが正方向に溢れる場合
				ret.tv_nsec -= ONE_SEC_IN_NANO;	// ナノ秒カウンタから1秒分引いて，
				ret.tv_sec++;					// 代わりに秒カウンタを加算
			}
			return ret;
		}
		
		//! @brief timespec構造体用の時刻を減算する関数 a - b
		//! @param[in]	a	timespec時刻
		//! @param[in]	b	timespec時刻
		//! @return 減算結果
		static timespec timespec_sub(const timespec& a, const timespec& b){
			timespec ret;
			ret.tv_sec = a.tv_sec - b.tv_sec;		// aとbの秒カウンタを減算
			ret.tv_nsec = a.tv_nsec - b.tv_nsec;	// aとbのナノ秒カウンタを減算
			if(ret.tv_nsec < 0){
				// ナノ秒カウンタが負方向に溢れる場合
				ret.tv_nsec += ONE_SEC_IN_NANO;	// ナノ秒カウンタに1秒分加えて，
				ret.tv_sec--;					// 代わりに秒カウンタを減算
			}
			return ret;
		}
		
		//! @brief timespec構造体用の時刻を比較する関数 a <= b
		//! @param[in]	a	timespec時刻
		//! @param[in]	b	timespec時刻
		//! @return 比較結果
		static bool timespec_lessthaneq(const timespec& a, const timespec& b){
			if(a.tv_sec != b.tv_sec){
				if(a.tv_sec < b.tv_sec){
					return true;
				}else{
					return false;
				}
			}else{
				if(a.tv_nsec <= b.tv_nsec){
					return true;
				}else{
					return false;
				}
			}
		}
		
		//! @brief timespec構造体用の時刻をクリアする関数 a = 0
		//! @param[in]	a	timespec時刻
		static void timespec_clear(timespec& a){
			a.tv_sec  = 0;	// 秒カウンタをクリア
			a.tv_nsec = 0;	// ナノ秒カウンタをクリア
		}
		
		//! @brief ナノ秒からtimespec構造体へ変換する関数
		//! @param[in]	時刻 [ns]
		//! @return		時刻(timespec構造体)
		static timespec nsec_to_timespec(const unsigned long time){
			timespec ret;
			ret.tv_sec = time/ONE_SEC_IN_NANO;
			ret.tv_nsec = time - ret.tv_sec*ONE_SEC_IN_NANO;
			return ret;
		}
		
		//! @brief timespec構造体から秒へ変換する関数
		//! @param[in]	時刻(timespec構造体)
		//! @return		時刻 [s]
		static double timespec_to_sec(const timespec& time){
			return (double)time.tv_sec + (double)time.tv_nsec*1e-9;	// [s] 時刻に変換して返す
		}
		
		//! @brief カーネルパラメータをリアルタイム用に設定する関数
		void SetKernelParameters(void){
			// WSLの場合
			if(IsInWSL() == true){
				EventLog("Setting kernel parameters for WSL");
				return;	// WSLのときは何もせずに終了
			}
			
			// x86_64系の場合
			#ifdef __x86_64__
				EventLog("Setting kernel parameters for x86_64");
				LinuxCommander::Execute("/bin/dmesg -n 1");	// カーネルメッセージの表示を抑制
				
				// カーネルパラメータの設定(CFSを無効にする場合)
				if constexpr( static_cast<uint8_t>(SFK) & static_cast<uint8_t>(SFkernelparam::CFS_DISABLED) ){
					// 下記が無いと1秒毎にスパイク状の遅延が生じる
					LinuxCommander::Execute("/bin/echo -1 > /proc/sys/kernel/sched_rt_runtime_us");			// CFSを無効
					LinuxCommander::Execute("/bin/echo 2147483647 > /proc/sys/kernel/sched_rt_period_us");	// リアルタイムタスク割り当て時間を最大化
					EventLog("Kernel param: CFS Disabled");
				}
				
				// カーネルパラメータの設定(PREEMPT_DYNAMICの場合にFULLモードに設定する場合)
				if constexpr( static_cast<uint8_t>(SFK) & static_cast<uint8_t>(SFkernelparam::PREEMPT_DYNFULL) ){
					LinuxCommander::Execute("/bin/echo \"full\" > /sys/kernel/debug/sched/preempt");
					EventLog("Kernel param: Set to \"full\" mode on PREEMPT_DYNAMIC");
				}
				
				// リアルタイムアルゴリズムのタイプによる設定
				if constexpr(SFA == SFalgorithm::WITHOUT_ZEROSLEEP){
					LinuxCommander::Execute("/bin/echo 0 > /proc/sys/kernel/watchdog");					// 「BUG: soft lockup」警告防止
					EventLog("Kernel param: Watchdog Disabled");
				}
			#endif
			
			// ARM系の場合
			#ifdef __ARM_ARCH
				EventLog("Setting kernel parameters for ARM");
				// パラメータなし
			#endif
			
			// 下記は実験的なカーネルパラメータ(様子見中)
			//LinuxCommander::Execute("/bin/echo 0 > /proc/sys/kernel/timer_migration");			// タイマの移行を無効
			//LinuxCommander::Execute("/bin/echo 0 > /proc/sys/kernel/sched_nr_migrate");			// プロセッサ間を移動できるタスク数をゼロにする
		}
		
		//! @brief カーネルパラメータを元に戻す関数
		void ResetKernelParameters(void){
			// WSLの場合
			if(IsInWSL() == true){
				return;	// WSLのときは何もせずに終了
			}
			
			// x86_64系の場合
			#ifdef __x86_64__
				// カーネルパラメータをデフォルト値に戻す(CFSを無効にした場合)
				if constexpr( static_cast<uint8_t>(SFK) & static_cast<uint8_t>(SFkernelparam::CFS_DISABLED) ){
					// 下記のように戻しておいた方が安定性の観点から無難
					LinuxCommander::Execute("/bin/echo 950000 > /proc/sys/kernel/sched_rt_runtime_us");	// CFSを有効，もとに戻す
					LinuxCommander::Execute("/bin/echo 1000000 > /proc/sys/kernel/sched_rt_period_us");	// リアルタイムタスク割り当て時間を元に戻す
					EventLog("Kernel param: CFS Enabled");
				}
				
				// カーネルパラメータをデフォルト値に戻す(PREEMPT_DYNAMICの場合にFULLモードに設定した場合)
				if constexpr( static_cast<uint8_t>(SFK) & static_cast<uint8_t>(SFkernelparam::PREEMPT_DYNFULL) ){
					LinuxCommander::Execute("/bin/echo \"voluntary\" > /sys/kernel/debug/sched/preempt");
					EventLog("Kernel param: Set to \"voluntary\" mode on PREEMPT_DYNAMIC");
				}
				
				LinuxCommander::Execute("/bin/dmesg -n 3");	// カーネルメッセージの表示をデフォルト値に戻す
				EventLog("Kernel param: Returned to default settings.");
			#endif
			
			// ARM系の場合
			#ifdef __ARM_ARCH
				// パラメータなし
			#endif
			
			// 下記は実験的なカーネルパラメータ(様子見中)
			//LinuxCommander::Execute("/bin/echo 1 > /proc/sys/kernel/watchdog");					// 「BUG: soft lockup」警告防止の解除(ウオッチドックは有効に戻せないらしい)
			//LinuxCommander::Execute("/bin/echo 1 > /proc/sys/kernel/timer_migration");			// タイマの移行を有効
			//LinuxCommander::Execute("/bin/echo 32 > /proc/sys/kernel/sched_nr_migrate");			// プロセッサ間を移動できるタスク数を既定値に戻す
		}
		
		//! @brief Windows Subsystem for Linux 内で動いてるのかチェックする関数
		//! @return	true = WSLで動いている，false = WSLではない
		static bool IsInWSL(void){
			const std::string WSLpath = "/proc/sys/fs/binfmt_misc/WSLInterop";	// WSL内に存在するあるファイル
			std::ifstream ifs(WSLpath);
			return ifs.is_open();
		}
};
}

#endif

