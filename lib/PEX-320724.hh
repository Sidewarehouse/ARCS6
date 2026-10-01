//! @file PEX320724.hh
//! @brief PEX-320724入出力クラス
//!
//! Interface社製入出力クラスのための入出力機能を提供します。(暫定短納期実装)
//!
//! @date 2026/10/01
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2026 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef PEX_320724
#define PEX_320724

#include <cassert>
#include "fbiad.h"

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
//! @brief PEX-320724入出力クラス
//! @tparam 
//template <>
class PEX320724 {
	public:
		//! @brief コンストラクタ
		PEX320724() noexcept
			: AdSmplConfig({}), ulSmplNum(0), SmplData{0}
		{
			// ADC開始
			int ret = AdOpen(1);
			arcs_assert(ret == AD_ERROR_SUCCESS);

			// サンプリング条件の取得
			AdBmGetSamplingConfig( 1, &AdSmplConfig );		// デフォルト設定を読み込み
			printf("fSmplFreq = %f Hz\n", AdSmplConfig.fSmplFreq);	// [Hz] デフォルトサンプリング周波数を表示
			AdSmplConfig.ulChCount = 1;						// チャンネル数
			AdSmplConfig.SmplChReq[0].ulChNo = 1;			// チャンネル番号1の設定
			AdSmplConfig.SmplChReq[0].ulRange = AD_5V;		// 入力電圧の範囲
			//AdSmplConfig.SmplChReq[1].ulChNo = 2;			// チャンネル番号2の設定
			//AdSmplConfig.SmplChReq[1].ulRange = AD_5V;		// 入力電圧の範囲
			AdSmplConfig.ulSingleDiff = AD_INPUT_SINGLE;	// シングルエンド入力に設定
			AdSmplConfig.ulSmplNum = 1024;					// [-] サンプリングするデータ数
			AdSmplConfig.ulSmplEventNum = 0;				// [-] 通知サンプリング件数
			AdSmplConfig.ulSmplRepeat = 1;					// [-] サンプリング繰り返し回数
			AdSmplConfig.ulBufferMode = AD_OVERWRITE;		// サンプリングデータ保存方法を上書き保存モードに設定
			AdSmplConfig.fSmplFreq = 1536000.0;				// [Hz] サンプリング周波数を1.536MHzに設定
			//AdSmplConfig.fSmplFreq = 48000.0;				// [Hz] サンプリング周波数を48kHzに設定
			AdSmplConfig.fScanFreq = 0;						// [Hz] チャンネルマルチプレクサ切替周波数
			AdSmplConfig.ulStartMode = AD_FREERUN;			// 開始トリガなしに設定
			AdSmplConfig.ulStopMode = AD_FREERUN;			// 終了トリガなしに設定
			AdSmplConfig.ulPreTrigDelay = 0;				// プリトリガディレイ
			AdSmplConfig.ulPostTrigDelay = 0;				// ポストトリガディレイ
			AdSmplConfig.ulATrgMode = AD_DISABLE;			// アナログトリガレベル、パルス出力無し
			AdSmplConfig.ulATrgPulse = AD_LOW_PULSE;		// アナログトリガパルス極性
			AdSmplConfig.ulStartTrigEdge = AD_DOWN_EDGE;	// 外部開始トリガ極性
			AdSmplConfig.ulStopTrigEdge = AD_DOWN_EDGE;		// 外部終了トリガ極性
			AdSmplConfig.ulTrigDI = 1;						// 汎用デジタル入力マスク設定
			AdSmplConfig.ulEClkEdge = AD_DOWN_EDGE;			// 外部クロック極性
			AdSmplConfig.ulFastMode = AD_NORMAL_MODE;		// 倍速モード使用しない
			AdSmplConfig.ulStatusMode = AD_NO_STATUS;		// サンプリングステータス
			AdSmplConfig.ulErrCtrl = AD_FREERUN;			// サンプリングエラー発生時の処理は無しに設定

			// サンプリング条件の設定
			ret = AdBmSetSamplingConfig(1, &AdSmplConfig);
			arcs_assert(ret != (int)AD_ERROR_NOW_SAMPLING);
			arcs_assert(ret != (int)AD_ERROR_INVALID_PARAMETER);
			arcs_assert(ret == AD_ERROR_SUCCESS);
			
			// 変換前準備
			AdSetFilter(1, AD_DF_8);		// オーバサンプリンレート
			AdClearSamplingData(1);			// サンプリングバッファクリア
		}
		
		//! @brief デストラクタ
		~PEX320724() noexcept {
			// ADC終了
			AdClose(1);
		}

		//! @brief 動作確認試験用関数
		void Test(void){
			// 変換開始＆ブロッキング動作＆データ取得
			AdStartSampling(1, FLAG_SYNC);
			unsigned long ulSmplNum = 1024;
			AdGetSamplingData(1, &SmplData[0], &ulSmplNum);	// データの取得

			// サンプリング動作状態の取得
			unsigned long ulAdSmplStatus = 0;
			unsigned long ulAdSmplCount = 0;
			unsigned long ulAdAvailCount = 0;
			AdGetStatus(1, &ulAdSmplStatus, &ulAdSmplCount, &ulAdAvailCount);
			printf("Status: 0x%08lx, N = %ld\n", ulAdSmplStatus, ulAdSmplCount);

			// 変換データの表示
			for(size_t i = 0; i < 1024; ++i) printf("%5ld: %8x, %8.5f V\n", i, SmplData[i], static_cast<double>(SmplData[i])*10.0/16777216.0 - 5.0);			
		}
		
		//! @brief 電圧を取得する関数
		//! @param[out]	V1	[V] CH1の電圧
		void GetVoltage(double& V1){
			// 変換開始＆ブロッキング動作＆データ取得
			AdStartSampling(1, FLAG_SYNC);
			unsigned long ulSmplNum = 1024;
			AdGetSamplingData(1, &SmplData[0], &ulSmplNum);	// データの取得

			// 1制御周期中の平均値を計算（オーバサンプリングと同等の動作）
			double Vave1 = 0;
			for(size_t i = 0; i < Nsmpl; ++i){
				Vave1 += static_cast<double>(SmplData[i])*10.0/16777216.0 - 5.0;
			}
			V1 = Vave1/static_cast<double>(Nsmpl);
		}

	private:
		PEX320724(const PEX320724&) = delete;					//!< コピーコンストラクタ使用禁止
		const PEX320724& operator=(const PEX320724&) = delete;	//!< コピー代入演算子使用禁止
		static constexpr size_t Nsmpl = 1024;	// [-] 連続AD変換するデータ数

		ADBMSMPLREQ AdSmplConfig;		// ADC設定用構造体
		unsigned long ulSmplNum;		// [-] サンプリングするデータ数
		uint32_t SmplData[1024];		// 連続サンプリングデータ
};
}

#endif

