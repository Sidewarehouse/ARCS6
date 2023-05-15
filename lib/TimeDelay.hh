//! @file TimeDelay.hh
//! @brief 遅延器クラス
//!
//! 遅延器 G(s)=e^(-s*T*n) -> G(z)=z^(-n)
//! (但し T は制御周期で，制御周期単位でしか遅延時間を設定できないことに注意)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TIMEDELAY
#define TIMEDELAY

namespace ARCS {	// ARCS名前空間
	//! @brief 遅延器クラス
	class TimeDelay {
		public:
			explicit TimeDelay(const long MaxDelay);// コンストラクタ MaxDelay；最大遅延時間
			~TimeDelay();							// デストラクタ
			double GetSignal(const double u);		// 出力信号の取得 u；入力信号
			void SetDelayTime(const long DelayTime);//遅延時間の設定 DelayTime；遅延時間 (最大遅延時間を越えないこと)
			void ClearDelayMemory(void);			// 遅延メモリのゼロクリア
			
		private:
			TimeDelay(const TimeDelay&);					// コピーコンストラクタ使用禁止
			const TimeDelay& operator=(const TimeDelay&);	// 代入演算子使用禁止
			long dmem_max;	// 最大遅延時間
			long num;		// 遅延時間
			long Wcount;	// 書き込みカウンタ
			long Rcount;	// 読み出しカウンタ
			double* dmem;	// 遅延メモリ
	};
}

#endif



