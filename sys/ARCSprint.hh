//! @file ARCSprint.hh
//! @brief ARCS デバッグプリントクラス
//!
//! ARCSデバッグ用のprintクラス
//! 好きな場所からARCS画面のデバッグ表示ができる。
//! ARCS6からリアルタイム空間でも使用可能になりました。
//!
//! @date 2022/02/20
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSPRINT
#define ARCSPRINT

#include <string>
#include "ARCScommon.hh"
#include "Matrix.hh"

// 関数呼び出し用マクロ
#define DebugPrint(a) (ARCSprint::DebugPrint_from_macro(a))						//!< 文字列用デバッグプリントマクロ  str : 表示する文字列
#define DebugPrintVar(a) (ARCSprint::DebugPrintVar_from_macro(a,#a))			//!< 実数用デバッグプリントマクロ  a : 表示する変数
#define DebugPrintVarFmt(a,b) (ARCSprint::DebugPrintVarFmt_from_macro(a,b,#a))	//!< 実数用デバッグプリントマクロ  a : 表示する変数, b : 表示形式 (%1.3e とか %5.3f とか printfと同じ)
#define DebugPrintMat(a) (ARCSprint::DebugPrintMat_from_macro(a,#a))			//!< 実数用デバッグプリントマクロ  a : 表示する変数
#define DebugIndic(a,b) (ARCSprint::DebugIndicator_from_macro(a,b))				//!< デバッグインジケータマクロ a : 表示する変数

// 前方宣言
namespace ARCS {
	class ARCSscreen;
}

namespace ARCS {
//! @brief ARCS printライブラリ
class ARCSprint {
	public:
		ARCSprint(void);	//!< コンストラクタ
		~ARCSprint();		//!< デストラクタ
		void SetScreenPtr(ARCSscreen* ScrPtr);	//!< ARCS画面ポインタの設定
		
		// 下記の関数はマクロから呼ばれることを想定
		static void DebugPrint_from_macro(const std::string& str);								//!< 文字列用デバッグプリント  str : 表示する文字列
		static void DebugPrintVar_from_macro(const double u, const std::string& varname);		//!< 変数用デバッグプリント  u : 表示する変数, varname : 変数名
		static void DebugPrintVarFmt_from_macro(const double u, const std::string& format, const std::string& varname);	//!< 変数用デバッグプリント  u : 表示する変数, format : 表示形式 (%1.3e とか %5.3f とか printfと同じ), varname : 変数名
		static void DebugIndicator_from_macro(const uint64_t Value, const unsigned int VarNo);	//!< デバッグインジケータへの値の設定
		
		//! @brief 行列用デバッグプリント
		//! @param[in] u 表示する行列
		//! @param[in] varname 変数名
		template <size_t NN, size_t MM>
		static void DebugPrintMat_from_macro(const Matrix<NN,MM>& u, const std::string& varname){
			size_t N = u.GetWidthLength();			// 行列の幅を取得
			size_t M = u.GetHeightLength();			// 行列の高さを取得
			std::string strbuff;					// 数値→文字列変換バッファ
			DebugPrint_from_macro(varname + " = ");	// 変数名を表示
			
			for(size_t j = 1; j <= M; ++j){
				strbuff += "[ ";
				for(size_t i = 1; i <= N; ++i){
					strbuff += ARCScommon::DoubleToString(u.GetElem(j,i), "% g");	// 行列要素を文字列へ変換
					if(i<N)strbuff += "  ";	// 要素間は空白で分離
				}
				strbuff += " ]";
				DebugPrint_from_macro(strbuff);	// ここで行を表示
				strbuff = "";			// バッファクリア
			}
		}
		
	private:
		ARCSprint(const ARCSprint&) = delete;					// コピーコンストラクタ使用禁止
		const ARCSprint& operator=(const ARCSprint&) = delete;	// 代入演算子使用禁止
		
		static ARCSscreen* ARCSscreenPtr;	//!< ARCS画面ポインタ
};
}

#endif

