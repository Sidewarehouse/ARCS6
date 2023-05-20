//! @file TwoInertiaStateFeedback.cc
//! @brief 2慣性共振系の状態フィードバックレギュレータ
//!
//! 2慣性共振系に対して，状態オブザーバの推定値を用いて状態フィードバックを掛けて任意の極配置にセットするレギュレータ。
//! (MATLABでいうところの「reg」のようなもの)
//!
//! @date 2020/04/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TwoInertiaStateFeedback.hh"

// テンプレートクラスのため，実体もヘッダ側に実装。

