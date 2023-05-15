//! @file Signum.cc
//! @brief 符号関数
//! 
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "ARCSassert.hh"
#include "Signum.hh"

//! @brief 符号関数
//! @param[in]	input	入力信号
//! @return 入力が正のとき1，0のとき0，負のとき-1
double ARCS::sgn(double input){
	if(0 < input)  return  1;	// 正のとき
	if(input == 0) return  0;	// 0のとき
	if(input < 0)  return -1;	// 負のとき
	arcs_assert(false);	// 例外
	return 0;
}


