//! @file ARCS.cc
//! @brief Advanced Robot Control System V6, ARCS6
//! @date 2024/06/26
//! @author Yokokura, Yuki
//!
//! @par コーディングの信条
//! - 低ジッタで低レイテンシのリアルタイム性をなによりも重視したコードを書くこと。
//! - 一方で、実行速度の過度な追求よりも可読性を重視したコードを書くこと。
//! - C++17をベースにしているが，古典的なC99の記述も許容しても良い。
//! - コードに文句があるなら，口を出した後に，手も動かすこと。
//! - 一貫性，対称性，双対性を考慮せよ。
//! - コンパイラの最適化に頼って良い。
//! - いきなり100点を目指すな。60点で動いてから改良せよ。
//!
//! @par コーディング既約
//! - システムコードを改変した場合は ./sys/ARCSparams.hh の ARCSparams::ARCS_REVISION を更新すること。
//! - 一時的にコメントアウトしたコードは最終的に消去すること。
//! - using namespace std; は例外なく使用禁止。
//! - 定数値は #define ではなく constexpr を使用すること。マクロは使って良い。
//! - グローバル変数の使用は極力回避すること。
//! - goto文の使用は極力回避すること。
//! - コメント文は「口語的表現」でなるべく付けて，コードとの相違がないようにすること。
//! - プロトタイプ宣言の引数変数は省略せずに書くこと。
//! - Efficient C++, More Efficient C++ にできるだけ準拠すること。
//! - コードを書いたらCppcheck等々により静的解析をすること。
//! - 物理定数，物理変数はすべてSI単位系で記述すること。
//!
//! @par エラー処理について
//! - 本プログラムは全体的に亘ってエラー処理・例外処理は行っていない。動かないときは動かない。
//!   ただし，arcs_assert関数を徹底的に埋め込んで常に動作チェックすること。(開発者:=ユーザというコンセプト)
//!
//! Copyright (C) 2011-2024 Yokokura, Yuki
//! MIT License. For details, see the LICENSE file.

#include <unistd.h>
#include "ARCScommon.hh"
#include "ARCSscreen.hh"
#include "ARCSthread.hh"
#include "ARCSeventlog.hh"
#include "ARCSassert.hh"
#include "ARCSprint.hh"
#include "ARCSparams.hh"
#include "ARCSscrparams.hh"
#include "ARCSgraphics.hh"

using namespace ARCS;

//! @brief ARCS エントリポイント
int main(void){
	// ここがすべての始まり(エントリポイント)
	// 初期設定
	ARCSeventlog ARCSlog;	// イベントログの初期化
	ARCSassert ARCSast;		// ARCS用assertの初期化
	ARCSprint ARCSprt;		// デバッグプリントの初期化
	PassedLog();			// イベントログにココを通過したことを記録
	
	// main関数のCPUコアとポリシーの設定
	ARCScommon::SetCPUandPolicy(pthread_self(), ARCSparams::ARCS_CPU_MAIN, ARCSparams::ARCS_POL_MAIN, ARCSparams::ARCS_PRIO_MAIN);
	
	ARCSscrparams ARCSscrpar;	// ARCS画面パラメータの生成
	ARCSgraphics ARCSgrph;		// ARCSグラフプロットの生成
	ARCSscreen ARCSscr(ARCSlog, ARCSast, ARCSprt, ARCSscrpar, ARCSgrph);// ARCS画面初期化＆初期画面描画
	ARCSthread ARCSthd(ARCSast, ARCSscrpar, ARCSgrph);					// リアルタイムスレッドの生成
	
	// 開始 or 終了指令入力の待機
	if(ARCSscr.WaitStartOrExit() == ARCSscreen::PHAS_EXIT){		// 「EXIT」が押された場合は，
		PassedLog();		// イベントログにココを通過したことを記録
		return EXIT_SUCCESS;// ARCS全終了
	}
	
	// EXIT系ボタンが押されるまでループ
	while(1){
		// 制御ここから
		PassedLog();		// イベントログにココを通過したことを記録
		ARCSthd.Start();	// リアルタイムスレッドを開始
		ARCSscr.WaitStop();	// 停止するまで待機(ブロッキング)
		ARCSthd.Stop();		// リアルタイムスレッドを停止
		PassedLog();		// イベントログにココを通過したことを記録
		// 制御ここまで
		
		// 再開始 or 終了指令入力の待機
		enum ARCSscreen::PhaseStatus PhaseRef = ARCSscr.WaitRestartOrDiscOrSaveExit();
		
		// 指令に従って挙動を変更
		if(PhaseRef == ARCSscreen::PHAS_RESTART){			// 「RESTART」が押された場合は，
			ARCSthd.Reset();								// リアルタイムスレッドをリセットしてループ続行
		}else if(PhaseRef == ARCSscreen::PHAS_SAVEEXIT){	// 「SAVE and EXIT」が押された場合は，
			ARCSthd.SaveDataFiles();						// 測定データを保存して，
			break;											// ループを抜ける
		}else{
			break;	// 「DISCARD and EXIT」が押された場合は何もせずループを抜ける
		}
	}
	
	PassedLog();		// イベントログにココを通過したことを記録
	return EXIT_SUCCESS;// ARCS全終了
	// ここがすべての終わり
}

