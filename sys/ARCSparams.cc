//! @file ARCSparams.hh
//! @brief ARCSシステムコード共通パラメータ設定静的クラス
//!
//! ARCSシステムコード内で共通に使用するパラメータ設定のための静的関数クラス
//!
//! @date 2024/05/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include "ARCSparams.hh"

using namespace ARCS;

// ARCS改訂番号(ARCS本体側システムコード改変時にちゃんと変えること)
const std::string ARCSparams::ARCS_REVISION("AR6-REV.24050214");	//!< (16文字以内)

// イベントログの設定
const std::string ARCSparams::EVENTLOG_NAME("EventLog.txt");		//!< イベントログファイル名

