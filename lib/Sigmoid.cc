//! @file Sigmoid.cc
//! @brief シグモイド関数
//! y = 1/(1 + exp(-a*u))
//! @date 2017/06/09
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2017 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "Sigmoid.hh"
#include <cmath>

//! @brief シグモイド関数
//! @param[in] u 入力信号
//! @param[in] a ゲイン
//! @return 出力信号
double ARCS::Sigmoid(double u, double a){
	return 1.0/(1.0 + exp(-a*u));
}


