//! @file SquareWave.cc
//! @brief 方形波発生器
//!
//! @date 2022/03/18
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cmath>
#include "SquareWave.hh"

//! @brief 方形波を出力する関数
//! @param[in]	freq	[Hz] 周波数
//! @param[in]	phase	[rad]位相
//! @param[in]	time	[s]  時刻
//! @return	方形波
double ARCS::SquareWave(const double freq, const double phase, const double time){
	double y;
	double r = sin(2.0*M_PI*freq*time + phase);
	if(0 < r){
		y =  1;
	}else{
		y = -1;
	}
	return y;
}

//! @brief 方形波を出力する関数（開始時刻指定版）
//! @param[in]	freq		[Hz] 周波数
//! @param[in]	phase		[rad]位相
//! @param[in]	time		[s]  時刻
//! @param[in]	starttime	[s]  開始時刻
//! @return	方形波
double ARCS::SquareWave(const double freq, const double phase, const double time, const double starttime){
	if(time < starttime) return 0;						// 開始時刻より前のときはゼロ出力
	return SquareWave(freq, phase, time - starttime);	// 開始時刻以降は方形波出力
}
