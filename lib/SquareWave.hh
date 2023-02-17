//! @file SquareWave.cc
//! @brief 方形波発生器
//!
//! @date 2022/01/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SQUAREWAVE
#define SQUAREWAVE

#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
double SquareWave(const double freq, const double phase, const double time);	//!< 方形波を出力する関数
double SquareWave(const double freq, const double phase, const double time, const double starttime);	//!< 方形波を出力する関数（開始時刻指定版）

//! @brief 方形波を出力する関数(縦ベクトル版)
//! @tparam	M	縦ベクトルの高さ
//! @param[in]	freq	[Hz] 周波数
//! @param[in]	phase	[rad]位相
//! @param[in]	time	[s]  時刻
//! @param[out]	output	方形波ベクトル
template <size_t M>
void SquareWave(const double freq, const double phase, const double time, Matrix<1,M>& output){
	double r = sin(2.0*M_PI*freq*time + phase);
	if(0 < r){
		output =  Matrix<1,M>::ones();
	}else{
		output = -Matrix<1,M>::ones();
	}
}

//! @brief 方形波を出力する関数(縦ベクトル版，ベクトル返し版)
//! @tparam	M	縦ベクトルの高さ
//! @param[in]	freq	[Hz] 周波数
//! @param[in]	phase	[rad]位相
//! @param[in]	time	[s]  時刻
//! @return	方形波ベクトル
template <size_t M>
Matrix<1,M> SquareWaveVec(const double freq, const double phase, const double time){
	Matrix<1,M> ret;
	SquareWave(freq, phase, time, ret);
	return ret;
}

}

#endif
