//! @file ARCSthread.hh
//! @brief ARCSリアルタイムスレッド管理クラス
//!
//! リアルタイムスレッドの生成、開始、停止、破棄などの管理をします。
//!
//! @date 2020/04/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ARCSTHREAD
#define ARCSTHREAD

#include <pthread.h>
#include <memory>
#include <functional>
#include "ControlFunctions.hh"
#include "SFthread.hh"
#include "DataMemory.hh"

// 前方宣言
namespace ARCS {
	class ARCSassert;
	class ScreenParams;
	class GraphPlot;
}

namespace ARCS {	// ARCS名前空間
	//! @brief ARCSリアルタイムスレッド管理クラス
	class ARCSthread {
		public:
			ARCSthread(ARCSassert& ARCSast, ScreenParams& SP, GraphPlot& GP);	//!< コンストラクタ
			~ARCSthread();					//!< デストラクタ
			void Start(void);				//!< スレッドを開始する関数
			void Stop(void);				//!< スレッドを停止する関数
			void Reset(void);				//!< スレッドをリセットする関数
			void SaveDataFiles(void);		//!< 測定データを保存する関数
			
		private:
			ARCSthread(const ARCSthread&) = delete;					//!< コピーコンストラクタ使用禁止
			const ARCSthread& operator=(const ARCSthread&) = delete;//!< 代入演算子使用禁止
			
			ARCSassert& ARCSast;	//!< ARCSアサートへの参照
			ScreenParams& ScrPara;	//!< 画面パラメータへの参照
			GraphPlot& Graph;		//!< グラフプロットへの参照
			DataMemory ExpDatMem;	//!< 実験データ保存メモリ
			
			ControlFunctions CtrlFuncs;								//!< 制御用周期実行関数群
			std::array<std::function<bool(double,double,double)>, ConstParams::THREAD_MAX> CtrlFuncObj;			//!< 制御用周期実行関数の関数オブジェクト配列
			std::array<std::unique_ptr<SFthread<ConstParams::THREAD_TYPE>>, ConstParams::THREAD_MAX> RTthreads;	//!< リアルタイムマルチスレッドへのスマートポインタ配列
			
			//! @brief スレッド状態の定義
			enum InfoThreadState {
				ITS_IDLE,	//!< アイドル状態
				ITS_START,	//!< 開始状態
				ITS_DSTRCT	//!< 破棄
			};
			InfoThreadState InfoState;	//!< 情報取得スレッドの状態
			pthread_mutex_t	InfoMutex;	//!< 情報取得スレッド同期用Mutex
			pthread_cond_t InfoCond;	//!< 情報取得スレッド同期用条件
			pthread_t InfoGetThreadID;						//!< 情報取得スレッドの識別子
			static void InfoGetThread(ARCSthread* const p);	//!< 情報取得スレッド
	};
}

#endif

