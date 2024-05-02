//! @file ARCSparams.hh
//! @brief ARCSシステムコード共通パラメータ設定静的クラス
//!
//! ARCSシステムコード内で共通に使用するパラメータ設定のための静的関数クラス
//!
//! @date 2024/05/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSPARAMS
#define ARCSPARAMS

#include <pthread.h>
#include <string>

namespace ARCS {	// ARCS名前空間
	//! @brief ARCSシステムコード共通パラメータ設定静的クラス
	class ARCSparams {
		public:
			// ARCS改訂番号(ARCS本体側システムコード改変時にちゃんと変えること)
			static const std::string ARCS_REVISION;	//!< (16文字以内)
			
			// イベントログの設定
			static const std::string EVENTLOG_NAME;	//!< イベントログファイル名
			
			// ARCSシステムスレッドの設定
			static constexpr int ARCS_POL_CMDI = SCHED_RR;	//!< 指令入力スレッドのポリシー
			static constexpr int ARCS_POL_DISP = SCHED_RR;	//!< 表示スレッドのポリシー
			static constexpr int ARCS_POL_EMER = SCHED_RR;	//!< 緊急停止スレッドのポリシー
			static constexpr int ARCS_POL_GRPL = SCHED_RR;	//!< グラフ表示スレッドのポリシー
			static constexpr int ARCS_POL_INFO = SCHED_RR;	//!< 情報取得スレッドのポリシー
			static constexpr int ARCS_POL_MAIN = SCHED_RR;	//!< main関数のポリシー
			static constexpr int ARCS_PRIO_CMDI = 32;		//!< 指令入力スレッドの優先順位(SCHED_RRはFIFO+32にするのがPOSIX.1-2001での決まり)
			static constexpr int ARCS_PRIO_DISP = 33;		//!< 表示スレッドの優先順位
			static constexpr int ARCS_PRIO_EMER = 34;		//!< 緊急停止スレッドの優先順位
			static constexpr int ARCS_PRIO_GRPL = 35;		//!< グラフ表示スレッドの優先順位
			static constexpr int ARCS_PRIO_INFO = 36;		//!< 情報取得スレッドの優先順位
			static constexpr int ARCS_PRIO_MAIN = 37;		//!< main関数スレッドの優先順位
			static constexpr unsigned int  ARCS_CPU_CMDI = 0;		//!< 指令入力スレッドに割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned int  ARCS_CPU_DISP = 0;		//!< 表示スレッドに割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned int  ARCS_CPU_EMER = 0;		//!< 緊急停止スレッドに割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned int  ARCS_CPU_GRPL = 1;		//!< グラフ表示スレッドに割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned int  ARCS_CPU_INFO = 0;		//!< 情報取得スレッドに割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned int  ARCS_CPU_MAIN = 0;		//!< main関数に割り当てるCPUコア番号（実時間スレッドとは別にすること）
			static constexpr unsigned long ARCS_TIME_DISP = 33333;	//!< [us] 表示の更新時間（ここの時間は厳密ではない）
			static constexpr unsigned long ARCS_TIME_GRPL = 33333;	//!< [us] グラフ表示の更新時間（ここの時間は厳密ではない）
			static constexpr unsigned long ARCS_TIME_INFO = 33333;	//!< [us] 情報取得の更新時間（ここの時間は厳密ではない）
			
		private:
			ARCSparams() = delete;	//!< コンストラクタ
			~ARCSparams() = delete;	//!< デストラクタ
			ARCSparams(const ARCSparams&) = delete;					//!< コピーコンストラクタ使用禁止
			const ARCSparams& operator=(const ARCSparams&) = delete;//!< 代入演算子使用禁止
			
	};
}

#endif

