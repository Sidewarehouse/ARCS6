//! @file SFthread.cc
//! @brief SCHED_FIFOリアルタイムスレッドクラス(sleep使用不使用テンプレート可変版, 関数オブジェクト版)
//!
//! pthreadのSCHED_FIFOで実時間スレッドを生成＆管理＆破棄する。実際に計測された制御周期や計算消費時間も提供する。
//!
//! @date 2020/04/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "SFthread.hh"

using namespace ARCS;

// テンプレートクラスのため，実体もヘッダ側に実装。

