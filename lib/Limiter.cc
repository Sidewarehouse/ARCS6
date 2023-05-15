//! @file Limiter.cc
//! @brief リミッタ
//! 任意の数値で入力を制限して出力
//! @date 2019/08/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "Limiter.hh"

//! @brief リミッタ
//! @param[in]	input	入力信号
//! @param[in]	limit	制限値
//! @return 制限された出力信号
double ARCS::Limiter(const double input, const double limit){
	double ret = input;
	if(limit <  input) ret =  limit;		// 任意の数値で入力を制限して出力
	if(input < -limit) ret = -limit;
	return ret;
}

//! @brief リミッタ
//! @param[in]	input		入力信号
//! @param[in]	pos_limit	正側制限値
//! @param[in]	neg_limit	負側制限値
//! @return 制限された出力信号
double ARCS::Limiter(const double input, const double pos_limit, const double neg_limit){
	double ret = input;
	if(pos_limit < input) ret = pos_limit;	// 任意の数値で入力を制限して出力
	if(input < neg_limit) ret = neg_limit;	// 任意の数値で入力を制限して出力
	return ret;
}
