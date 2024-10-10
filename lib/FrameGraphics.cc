//! @file FrameGraphics.cc
//! @brief フレームグラフィックスクラスV2(新型テンプレート版)
//!
//! LinuxフレームバッファとPNG画像ファイルへのグラフィックス描画を行うクラス
//! (PNG画像描画のみなら Windows Subsystem for Linux でも実行可能)
//! 32bit色のみ対応
//!
//! 構造の概略：
//! 「フレームバッファ」「画面バッファ」「背景バッファ」の３つがある。
//! Draw系関数を呼ぶと画面バッファに書き込まれる。
//! RefreshFrame関数を呼ぶと画面バッファからフレームバッファに画像が転送され，ディスプレイに表示される。
//! 現在の画面バッファを背景バッファとして取っておいて，後で読み出すこともできる。
//! 画面バッファをPNG画像ファイルとして保存することも可能。
//! WSL上などフレームバッファが存在しないときはダミーのバッファを作成してやり過ごし，PNGファイルで出力する。
//!
//! @date 2024/10/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include "FrameGraphics.hh"

// テンプレートクラスのため，実体もヘッダ側に実装。

using namespace ARCS;
