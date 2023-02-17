//! @file TwoInertiaStateDistObsrv.hh
//! @brief 2慣性系状態外乱オブザーバクラス
//!
//! q軸電流とモータ側速度から，2慣性系の負荷側速度・ねじれ角・モータ側速度・負荷側ステップ外乱を推定する状態外乱オブザーバ
//! (普通の状態オブザーバと負荷側の0次外乱オブザーバを組み合わせたオブザーバ)
//!
//! @date 2020/04/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TwoInertiaStateDistObsrv.hh"

// テンプレートクラスのため，実体もヘッダ側に実装。
