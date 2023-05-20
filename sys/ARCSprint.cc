//! @file ARCSprint.hh
//! @brief ARCS デバッグプリントクラス
//!
//! ARCSデバッグ用のprintクラス
//! 好きな場所からARCS画面のデバッグ表示ができる。
//! ARCS6からリアルタイム空間でも使用可能になりました。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include "ARCSprint.hh"
#include "ARCSscreen.hh"

using namespace ARCS;

// 静的メンバ変数の実体
ARCSscreen* ARCSprint::ARCSscreenPtr = nullptr;	//!< ARCS画面ポインタ

//! @brief コンストラクタ
ARCSprint::ARCSprint(void){
	PassedLog();
}

//! @brief デストラクタ
ARCSprint::~ARCSprint(){
	PassedLog();
}

//! @brief ARCS画面ポインタの設定
//! @param[in]	ScrPtr	ARCS画面へのポインタ
void ARCSprint::SetScreenPtr(ARCSscreen* ScrPtr){
	ARCSscreenPtr = ScrPtr;
}

//! @brief 文字列用デバッグプリント
//! @param[in] str 表示する文字列
void ARCSprint::DebugPrint_from_macro(const std::string& str){
	// デバッグプリントを画面バッファに書き込む
	if(ARCSscreenPtr != nullptr){
		ARCSscreenPtr->WriteDebugPrintBuffer(str);
	}
}

//! @brief 変数用デバッグプリント
//! @param[in] u 表示する変数
//! @param[in] varname 変数名
void ARCSprint::DebugPrintVar_from_macro(const double u, const std::string& varname){
	DebugPrintVarFmt_from_macro(u, "% g", varname);
}

//! @brief 変数用デバッグプリント
//! @param[in] u 表示する変数
//! @param[in] format 表示形式 (%1.3e とか %5.3f とか printfと同じ)
//! @param[in] varname 変数名
void ARCSprint::DebugPrintVarFmt_from_macro(const double u, const std::string& format, const std::string& varname){
	std::string str;
	str = varname + " = " + ARCScommon::DoubleToString(u, format);	// 変数名と数値を文字列へ
	DebugPrint_from_macro(str);
}

//! @brief デバッグインジケータへの値の設定
//! @param[in]	Value	整数値
//! @param[in]	VarNo	変数番号
void ARCSprint::DebugIndicator_from_macro(const uint64_t Value, const unsigned int VarNo){
	// デバッグインジケータを画面バッファに書き込む
	if(ARCSscreenPtr != nullptr){
		ARCSscreenPtr->WriteDebugIndicator(Value, VarNo);
	}
}

