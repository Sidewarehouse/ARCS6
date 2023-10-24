//! @file FunctionBase.cc
//! @brief 関数群のベースコード
//!
//! クラスにまとめるまでもない関数群を追加する場合は，この関数ベースコードを基に作ってネ。
//!
//! @date 20XX/XX/XX
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <cassert>
#include "FunctionBase.hh"

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

using namespace ARCS;

//! @brief 関数のベース
//! @param[in]	u	入力
//! @return	出力
double ARCS::FunctionBase(double u){
	return u;
}
