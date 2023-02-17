//! @file TwoStairsWave.cc
//! @brief 2段階段波発生器
//!
//! @date 2022/11/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cmath>
#include "TwoStairsWave.hh"
#include "TriangleWave.hh"

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

//! @brief 2段階段波を出力する関数
//! @param[in]	freq	[Hz] 周波数
//! @param[in]	time	[s]  時刻
//! @return	方形波
double ARCS::TwoStairsWave(const double freq, const double time){
	double y;
	double r = TriangleWave(freq, time) + 1.0;	// 0～＋2の範囲の三角波
	if(0 <= r && r < 0.67){
		y = 0;
	}else if( 0.67 <= r && r < 1.33){
		y = 0.5;
	}else{
		y = 1;
	}
	return y;
}

//! @brief 2段階段波を出力する関数（開始時刻指定版）
//! @param[in]	freq		[Hz] 周波数
//! @param[in]	time		[s]  時刻
//! @param[in]	starttime	[s]  開始時刻
//! @return	方形波
double ARCS::TwoStairsWave(const double freq, const double time, const double starttime){
	if(time < starttime) return 0;					// 開始時刻より前のときはゼロ出力
	return TwoStairsWave(freq, time - starttime);	// 開始時刻以降は方形波出力
}
