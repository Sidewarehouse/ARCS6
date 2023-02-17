//! @file ARCScommon.hh
//! @brief ARCS共通静的関数クラス
//!
//! ARCS内で共通に使用する静的関数クラス
//!
//! @date 2020/04/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ARCSCOMMON
#define ARCSCOMMON

#include <string>

namespace ARCS {	// ARCS名前空間
	//! @brief ARCS共通静的関数クラス
	class ARCScommon {
		public:
			static void SetCPUandPolicy(const pthread_t FuncName, const unsigned int CPUnum, const int Policy, const int Priority); //!< スレッドを動作させるCPUとポリシーを設定する関数
			static std::string DoubleToString(double u, const std::string& format);		//!< 浮動小数点から文字列へ (std::to_stringだとフォーマット指定できないぽいので仕方なく実装)
			static std::string Uint64ToString(uint64_t u, const std::string& format);	//!< 整数値uint64_tから文字列へ
			static std::string GetNowTime(void);					//!< 現在時刻を取得する関数
			
		private:
			ARCScommon() = delete;	//!< コンストラクタ
			~ARCScommon() = delete;	//!< デストラクタ
			ARCScommon(const ARCScommon&) = delete;					//!< コピーコンストラクタ使用禁止
			const ARCScommon& operator=(const ARCScommon&) = delete;//!< 代入演算子使用禁止
			static constexpr unsigned int NUM_CHARBUFF = 32;		//!< 変換用文字列バッファ数
	};
}

#endif

