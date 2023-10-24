//! @file ARCSassert.hh
//! @brief ARCS用ASSERTクラス
//! @date 2023/10/23
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSASSERT
#define ARCSASSERT

#include <cassert>
#include <pthread.h>
#include <string>

// 関数呼び出し用マクロ(コンパイル時定数の場合assertを呼び出し、実行時の場合assert_from_macroを呼び出す)
#define arcs_assert(a) (__builtin_constant_p(a) ? assert(a) : ARCSassert::assert_from_macro(a,#a,__FILE__,__LINE__))	//!< ARCS用assertマクロ  a : assert条件

namespace ARCS {	// ARCS名前空間
	//! @brief ARCS用ASSERTクラス
	class ARCSassert {
		public:
			ARCSassert();	//!< コンストラクタ
			~ARCSassert();	//!< デストラクタ
			static void assert_from_macro(const bool assertion, const std::string& condition, const std::string& file, const int line); //!<ARCS用assert関数(マクロから呼ばれることを想定)
			bool IsRealtimeMode(void);		//!< リアルタイムモードかどうか返す関数
			void SetRealtimeMode(void);		//!< リアルタイムモードに切り替える関数
			void SetNonRealtimeMode(void);	//!< 非リアルタイムモードに切り替える関数
			bool IsEmergency(void);			//!< 緊急事態かどうか返す関数
			void DeclareEmergency(void);	//!< 緊急事態を宣言する関数
			void DoneEmergencyProc(void);	//!< 緊急停止処理が完了したことを知らせる関数
			bool IsDoneEmergencyProc(void);	//!< 緊急停止処理が完了したかどうか返す関数
			std::tuple< std::string, std::string, int >
			 GetAssertInfo(void);			//!< Assert条件，ファイル名，行番号を返す関数
			void WaitEmergency(void);		//!< 緊急停止が発生するまで待機する関数
			
		private:
			ARCSassert(const ARCSassert&) = delete;					//!< コピーコンストラクタ使用禁止
			const ARCSassert& operator=(const ARCSassert&) = delete;//!< 代入演算子使用禁止
			
			static const unsigned int WAIT_TIME = 1000;		//!< [us] ブロッキング待機の周期
			
			static pthread_mutex_t SyncMutex;		//!< Mutexロック
			static pthread_cond_t SyncCond;			//!< 待機条件
			
			static bool RealtimeModeFlag;			//!< リアルタイムモードフラグ
			static bool EmergencyStopFlag;			//!< 緊急停止フラグ  true = 緊急停止, false = 通常
			static bool EmergencyProcDoneFlag;		//!< 緊急停止処理完了フラグ true = 完了，false = 未完了
			static std::string EmergencyStopCond;	//!< 緊急停止時該当条件
			static std::string EmergencyStopFile;	//!< 緊急停止時ファイル名
			static int EmergencyStopLine;			//!< 緊急停止時行番号
	};
}

#endif

