//! @file Limiter.hh
//! @brief リミッタ
//!
//! 任意の数値で入力を制限して出力
//!
//! @date 2021/08/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef LIMITER
#define LIMITER

#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
	double Limiter(const double input, const double limit);	//!< リミッタ
	double Limiter(const double input, const double pos_limit, const double neg_limit);	//!< リミッタ
	
	//! @brief リミッタ(縦ベクトル版)
	//! @tparam	M	縦ベクトルの高さ
	//! @param[in,out]	inout	入出力ベクトル
	//! @param[in]		limit	リミット値ベクトル
	template<size_t M>
	void Limiter(Matrix<1,M>& inout, const Matrix<1,M>& limit){
		for(size_t i = 1; i <= M; ++i) inout[i] = Limiter(inout[i], limit[i]);
	}
}

#endif

