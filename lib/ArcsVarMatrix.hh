//! @file ArcsVarMatrix.hh
//! @brief ARCS-Variable-Matrix 可変行列演算クラス
//!
//!
//! @date 2025/12/31
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.
//

#ifndef ARCSVARMATRIX
#define ARCSVARMATRIX

#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdint>

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

#include "ArcsMatrix.hh"

// ARCS名前空間
namespace ARCS {


}

#endif
