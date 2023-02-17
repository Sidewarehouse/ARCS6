//! @file TwoInertiaStateObsrv.cc
//! @brief 2慣性系状態オブザーバクラス
//!
//! q軸電流とモータ側速度から，2慣性系の負荷側速度・ねじれ角・モータ側速度を推定する状態オブザーバ
//!
//! @date 2020/04/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TwoInertiaStateObsrv.hh"

// テンプレートクラスのため，実体もヘッダ側に実装。
