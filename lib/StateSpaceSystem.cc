//! @file StateSpaceSystem.cc
//! @brief 状態空間表現によるシステム
//!
//! 線形の状態空間モデルで表現されたシステムを保持，入力信号に対する出力信号を計算する。
//! (MATLABでいうところの「State Space」のようなもの)
//!
//! @date 2020/04/03
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "StateSpaceSystem.hh"

// テンプレートクラスのため，実体もヘッダ側に実装。
