//! @file StepWave.cc
//! @brief ステップ波形生成器
//!
//! ステップ波形を生成する。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include "StepWave.hh"

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

//! @brief 1回だけステップ波形を出力する関数
//! @param[in]	Tstp	ステップ時刻
//! @param[in]	t	現在時刻
//! @return	出力
double ARCS::StepWave(const double Tstp, const double t){
	double ret;
	if(t < Tstp){
		ret = 0;
	}else{
		ret = 1;
	}
	return ret;
}
