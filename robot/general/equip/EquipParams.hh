//! @file EquipParams.hh
//! @brief 実験装置用定数値格納用クラス
//!        ARCSに必要な実験装置に特有な定数値を格納します。
//! @date 2024/06/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef EQUIPPARAMS
#define EQUIPPARAMS

#include <cmath>
#include "ARCSparams.hh"
#include "SFthread.hh"

namespace ARCS {	// ARCS名前空間
//! @brief 実験装置用定数値格納用クラス
class EquipParams {
	public:
		// 画面サイズの設定 (モニタ解像度に合うように設定すること)
		// 1024×600(WSVGA) の場合に下記をアンコメントすること
		//static constexpr int SCR_VERTICAL_MAX = 36;		//!< [文字] 画面の最大高さ文字数
		//static constexpr int SCR_HORIZONTAL_MAX = 127;	//!< [文字] 画面の最大幅文字数
		// 1024×768(XGA) の場合に下記をアンコメントすること
		static constexpr int SCR_VERTICAL_MAX = 47;			//!< [文字] 画面の最大高さ文字数
		static constexpr int SCR_HORIZONTAL_MAX = 127;		//!< [文字] 画面の最大幅文字数
		// 1280×1024(SXGA) の場合に下記をアンコメントすること
		//static constexpr int SCR_VERTICAL_MAX = 63;		//!< [文字] 画面の最大高さ文字数
		//static constexpr int SCR_HORIZONTAL_MAX = 159;	//!< [文字] 画面の最大幅文字数
		// 1920×1080(Full HD) の場合に下記をアンコメントすること
		//static constexpr int SCR_VERTICAL_MAX = 66;		//!< [文字] 画面の最大高さ文字数
		//static constexpr int SCR_HORIZONTAL_MAX = 239;	//!< [文字] 画面の最大幅文字数
		// それ以外の解像度の場合は各自で値を探すこと
		
		// SCHED_FIFOリアルタイムスレッドの設定
		static constexpr SFalgorithm THREAD_TYPE = SFalgorithm::INSERT_ZEROSLEEP;	//!< リアルタイムアルゴリズムの選択
		// 上記を INSERT_ZEROSLEEP にすると安定性が増すがリアルタイム性は落ちる。遅い処理系の場合に推奨。
		// WITHOUT_ZEROSLEEP にするとリアルタイム性が向上するが，一時的に操作不能になる可能性が残る。高速な処理系の場合に選択可。
		// 下記はカーネルパラメータの設定
		// NO_SETTINGS もしくは CFS_DISABLED と PREEMPT_DYNFULL が併用可。詳細はSFthreadクラスのコメント欄を参照のこと。
		//static constexpr SFkernelparam THREAD_KP = SFkernelparam::CFS_DISABLED;		//!< CFSをリアルタイム用に設定
		static constexpr SFkernelparam THREAD_KP = static_cast<SFkernelparam>(
			static_cast<uint8_t>(SFkernelparam::CFS_DISABLED) | static_cast<uint8_t>(SFkernelparam::PREEMPT_DYNFULL)
		);	//!< CFSとPREEMPTの設定を併用する場合の例
		
		//! @brief 使用CPUコアの設定
		//! CPU0番コアはOSとARCSシステム、CPU1番コアはARCS描画系が使用しているので、2番目以上が望ましい
		static constexpr std::array<unsigned int, ARCSparams::THREAD_MAX> CPUCORE_NUMBER = {
				3,	// [-] 制御用周期実行関数1 (スレッド1) 使用するCPUコア番号
				2,	// [-] 制御用周期実行関数2 (スレッド2) 使用するCPUコア番号
				1,	// [-] 制御用周期実行関数3 (スレッド3) 使用するCPUコア番号
		};
		
		// 実験機アクチュエータの設定
		static constexpr unsigned int ACTUATOR_NUM = 1;		//!< [基] 実験装置のアクチュエータの総数
		
		//! @brief 実験機アクチュエータの種類の設定（リニアモータか回転モータかの設定）
		//! 下記が使用可能
		//! LINEAR_MOTOR = リニアモータ
		//!	ROTARY_MOTOR = 回転モータ
		static constexpr std::array<ARCSparams::ActType, ARCSparams::ACTUATOR_MAX> ACT_TYPE = {
			ARCSparams::ActType::ROTARY_MOTOR,	//  1番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  2番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  3番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  4番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  5番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  6番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  7番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  8番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	//  9番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 10番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 11番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 12番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 13番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 14番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 15番 アクチュエータ
			ARCSparams::ActType::ROTARY_MOTOR,	// 16番 アクチュエータ
		};
		
		//! @brief 実験機アクチュエータの指令単位の設定（電流なのか推力なのかトルクなのかの設定）
		//! 下記が使用可能
		//!	AMPERE = アンペア単位
		//!	NEWTON = ニュートン単位
		//!	NEWTON_METER = ニュートンメートル単位
		static constexpr std::array<ARCSparams::ActRefUnit, ARCSparams::ACTUATOR_MAX> ACT_REFUNIT = {
			ARCSparams::ActRefUnit::AMPERE,	//  1番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  2番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  3番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  4番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  5番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  6番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  7番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  8番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	//  9番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 10番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 11番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 12番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 13番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 14番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 15番 アクチュエータ
			ARCSparams::ActRefUnit::AMPERE,	// 16番 アクチュエータ
		};
		
		//! @brief トルク/推力定数の設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_FORCE_TORQUE_CONST = {
			1,	// [N/A]/[Nm/A]  1番 アクチュエータ
			1,	// [N/A]/[Nm/A]  2番 アクチュエータ
			1,	// [N/A]/[Nm/A]  3番 アクチュエータ
			1,	// [N/A]/[Nm/A]  4番 アクチュエータ
			1,	// [N/A]/[Nm/A]  5番 アクチュエータ
			1,	// [N/A]/[Nm/A]  6番 アクチュエータ
			1,	// [N/A]/[Nm/A]  7番 アクチュエータ
			1,	// [N/A]/[Nm/A]  8番 アクチュエータ
			1,	// [N/A]/[Nm/A]  9番 アクチュエータ
			1,	// [N/A]/[Nm/A] 10番 アクチュエータ
			1,	// [N/A]/[Nm/A] 11番 アクチュエータ
			1,	// [N/A]/[Nm/A] 12番 アクチュエータ
			1,	// [N/A]/[Nm/A] 13番 アクチュエータ
			1,	// [N/A]/[Nm/A] 14番 アクチュエータ
			1,	// [N/A]/[Nm/A] 15番 アクチュエータ
			1,	// [N/A]/[Nm/A] 16番 アクチュエータ
		};
		
		//! @brief 定格電流値の設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_RATED_CURRENT = {
			1,	// [A]  1番 アクチュエータ
			1,	// [A]  2番 アクチュエータ
			1,	// [A]  3番 アクチュエータ
			1,	// [A]  4番 アクチュエータ
			1,	// [A]  5番 アクチュエータ
			1,	// [A]  6番 アクチュエータ
			1,	// [A]  7番 アクチュエータ
			1,	// [A]  8番 アクチュエータ
			1,	// [A]  9番 アクチュエータ
			1,	// [A] 10番 アクチュエータ
			1,	// [A] 11番 アクチュエータ
			1,	// [A] 12番 アクチュエータ
			1,	// [A] 13番 アクチュエータ
			1,	// [A] 14番 アクチュエータ
			1,	// [A] 15番 アクチュエータ
			1,	// [A] 16番 アクチュエータ
		};
		
		//! @brief 瞬時最大許容電流値の設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_MAX_CURRENT = {
			3,	// [A]  1番 アクチュエータ
			3,	// [A]  2番 アクチュエータ
			3,	// [A]  3番 アクチュエータ
			3,	// [A]  4番 アクチュエータ
			3,	// [A]  5番 アクチュエータ
			3,	// [A]  6番 アクチュエータ
			3,	// [A]  7番 アクチュエータ
			3,	// [A]  8番 アクチュエータ
			3,	// [A]  9番 アクチュエータ
			3,	// [A] 10番 アクチュエータ
			3,	// [A] 11番 アクチュエータ
			3,	// [A] 12番 アクチュエータ
			3,	// [A] 13番 アクチュエータ
			3,	// [A] 14番 アクチュエータ
			3,	// [A] 15番 アクチュエータ
			3,	// [A] 16番 アクチュエータ
		};
		
		//! @brief 定格トルクの設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_RATED_TORQUE = {
			1,	// [Nm]  1番 アクチュエータ
			1,	// [Nm]  2番 アクチュエータ
			1,	// [Nm]  3番 アクチュエータ
			1,	// [Nm]  4番 アクチュエータ
			1,	// [Nm]  5番 アクチュエータ
			1,	// [Nm]  6番 アクチュエータ
			1,	// [Nm]  7番 アクチュエータ
			1,	// [Nm]  8番 アクチュエータ
			1,	// [Nm]  9番 アクチュエータ
			1,	// [Nm] 10番 アクチュエータ
			1,	// [Nm] 11番 アクチュエータ
			1,	// [Nm] 12番 アクチュエータ
			1,	// [Nm] 13番 アクチュエータ
			1,	// [Nm] 14番 アクチュエータ
			1,	// [Nm] 15番 アクチュエータ
			1,	// [Nm] 16番 アクチュエータ
		};
		
		//! @brief 瞬時最大トルクの設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_MAX_TORQUE = {
			3,	// [Nm]  1番 アクチュエータ
			3,	// [Nm]  2番 アクチュエータ
			3,	// [Nm]  3番 アクチュエータ
			3,	// [Nm]  4番 アクチュエータ
			3,	// [Nm]  5番 アクチュエータ
			3,	// [Nm]  6番 アクチュエータ
			3,	// [Nm]  7番 アクチュエータ
			3,	// [Nm]  8番 アクチュエータ
			3,	// [Nm]  9番 アクチュエータ
			3,	// [Nm] 10番 アクチュエータ
			3,	// [Nm] 11番 アクチュエータ
			3,	// [Nm] 12番 アクチュエータ
			3,	// [Nm] 13番 アクチュエータ
			3,	// [Nm] 14番 アクチュエータ
			3,	// [Nm] 15番 アクチュエータ
			3,	// [Nm] 16番 アクチュエータ
		};
		
		//! @brief 初期位置の設定
		static constexpr std::array<double, ARCSparams::ACTUATOR_MAX> ACT_INITPOS = {
			0,	// [rad]  1軸 アクチュエータ
			0,	// [rad]  2軸 アクチュエータ
			0,	// [rad]  3軸 アクチュエータ
			0,	// [rad]  4軸 アクチュエータ
			0,	// [rad]  5軸 アクチュエータ
			0,	// [rad]  6軸 アクチュエータ 
			0,	// [rad]  7番 アクチュエータ
			0,	// [rad]  8番 アクチュエータ
			0,	// [rad]  9番 アクチュエータ
			0,	// [rad] 10番 アクチュエータ
			0,	// [rad] 11番 アクチュエータ
			0,	// [rad] 12番 アクチュエータ
			0,	// [rad] 13番 アクチュエータ
			0,	// [rad] 14番 アクチュエータ
			0,	// [rad] 15番 アクチュエータ
			0,	// [rad] 16番 アクチュエータ
		};
		
	private:
		EquipParams() = delete;	//!< コンストラクタ使用禁止
		~EquipParams() = delete;//!< デストラクタ使用禁止
		EquipParams(const EquipParams&) = delete;					//!< コピーコンストラクタ使用禁止
		const EquipParams& operator=(const EquipParams&) = delete;	//!< 代入演算子使用禁止
};
}

#endif

