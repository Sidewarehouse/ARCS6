//! @file DeadBand.cc
//! @brief 不感帯
//! 
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "ARCSassert.hh"
#include "DeadBand.hh"

//! @brief 不感帯
//! @param[in]	u	入力信号
//! @param[in]	Upper	不感帯の上限
//! @param[in]	Lower	不感帯の下限
//! @return 不感帯通過後の信号
double ARCS::DeadBand(double u, double Upper, double Lower){
	if(Upper < u) return u - Upper;			// 不感帯よりも上のとき
	if(u < Lower) return u - Lower;			// 不感帯よりも下のとき
	return 0;	// 不感帯の中に入っているとき
}


