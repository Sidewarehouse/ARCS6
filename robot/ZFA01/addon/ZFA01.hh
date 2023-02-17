//! @file ZFA01.hh
//! @brief Zero Force Arm 01 - ZFA01クラス
//!
//! ZFA01ロボット制御のためのクラス
//!
//! @date 2022/03/29
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ZFA01_INCGUARD
#define ZFA01_INCGUARD

#include <cmath>
#include "Matrix.hh"
#include "TwoInertiaParamDef.hh"

namespace ARCS {	// ARCS名前空間
//! @brief Zero Force Arm 01 - ZFA01クラス
class ZFA01 {
	public:
		// ZFA01基本パラメータの定義
		static constexpr size_t N_AX = 3;				//!< [-] 軸数
		static constexpr double EPS_SINGULAR = 1e-4;	//!< [-] 特異点判定ゼロ閾値
		static constexpr double TORQUE_CONST1S = 0.320;	//!< [Nm/A] 1S軸 トルク定数 (ユニサーボ SVM-750-100 同定値)
		static constexpr double TORQUE_CONST2L = 0.372;	//!< [Nm/A] 2L軸 トルク定数 (ユニサーボ SVM-750-100 同定値)
		static constexpr double TORQUE_CONST3U = 0.393;	//!< [Nm/A] 3U軸 トルク定数 (ユニサーボ SVM-200-100 同定値)
		static constexpr double CURRENT_RATED1S = 7.62;	//!< [A]  1S軸 定格電流 (ユニサーボ SVM-750-100)
		static constexpr double CURRENT_RATED2L = 7.62;	//!< [A]  2L軸 定格電流 (ユニサーボ SVM-750-100)
		static constexpr double CURRENT_RATED3U = 2.25;	//!< [A]  3U軸 定格電流 (ユニサーボ SVM-200-100)
		static constexpr double CURRENT_MAX1S = 25.5;	//!< [A]  1S軸 瞬時最大電流 (ユニサーボ SVM-750-100)
		static constexpr double CURRENT_MAX2L = 25.5;	//!< [A]  2L軸 瞬時最大電流 (ユニサーボ SVM-750-100)
		static constexpr double CURRENT_MAX3U =  7.25;	//!< [A]  3U軸 瞬時最大電流 (ユニサーボ SVM-200-100)
		static constexpr double TORQUE_RATED1S = 174;	//!< [Nm]  1S軸 定格トルク (ユニサーボ SVM-750-100)
		static constexpr double TORQUE_RATED2L = 174;	//!< [Nm]  2L軸 定格トルク (ユニサーボ SVM-750-100)
		static constexpr double TORQUE_RATED3U =  45;	//!< [Nm]  3U軸 定格トルク (ユニサーボ SVM-200-100)
		static constexpr double TORQUE_MAX1S = 433;		//!< [Nm]  1S軸 瞬時最大トルク (ユニサーボ SVM-750-100)
		static constexpr double TORQUE_MAX2L = 433;		//!< [Nm]  2L軸 瞬時最大トルク (ユニサーボ SVM-750-100)
		static constexpr double TORQUE_MAX3U = 107;		//!< [Nm]  3U軸 瞬時最大トルク (ユニサーボ SVM-200-100)
		static constexpr double GEAR_RATIO1S = 100;		//!< [-] 1S軸 減速比 (ユニサーボ SVM-750-100)
		static constexpr double GEAR_RATIO2L = 100;		//!< [-] 2L軸 減速比 (ユニサーボ SVM-750-100)
		static constexpr double GEAR_RATIO3U = 100;		//!< [-] 3U軸 減速比 (ユニサーボ SVM-200-100)
		static constexpr double POSITION_INIT1S =  6.25931872328057e+00;				//!< [rad] 1S軸 初期位置
		static constexpr double POSITION_INIT2L =  5.97;								//!< [rad] 2L軸 初期位置
		static constexpr double POSITION_INIT3U =  6.90 + GEAR_RATIO3U*M_PI/180.0*135.0;//!< [rad] 3U軸 初期位置(初期角45°)
		static constexpr double POSITION_LIMPOS1 =   45.0*M_PI/180.0;	//!< [rad]	1S軸 正側可動範囲
		static constexpr double POSITION_LIMNEG1 = -180.0*M_PI/180.0;	//!< [rad]	1S軸 負側可動範囲
		static constexpr double POSITION_LIMPOS2 =   40.0*M_PI/180.0;	//!< [rad]	2L軸 正側可動範囲
		static constexpr double POSITION_LIMNEG2 = -100.0*M_PI/180.0;	//!< [rad]	2L軸 負側可動範囲
		static constexpr double POSITION_LIMPOS3 =    5.0*M_PI/180.0;	//!< [rad]	3U軸 正側可動範囲
		static constexpr double POSITION_LIMNEG3 = -137.0*M_PI/180.0;	//!< [rad]	3U軸 負側可動範囲
		static constexpr double INERTIA_ALL1S = 2.45e-4;		//!< [kgm^2] 1S軸 全慣性 Ja
		static constexpr double INERTIA_ALL2L = 3.28e-4;		//!< [kgm^2] 2L軸 全慣性 Ja
		static constexpr double INERTIA_ALL3U = 5.50e-5;		//!< [kgm^2] 3U軸 全慣性 Ja
		static constexpr double INERTIA_LOAD1S = 6.27e-2;		//!< [kgm^2] 1S軸 負荷側慣性
		static constexpr double INERTIA_LOAD2L = 1.20;			//!< [kgm^2] 2L軸 負荷側慣性
		static constexpr double INERTIA_LOAD3U = 0.249;			//!< [kgm^2] 3U軸 負荷側慣性
		static constexpr double INERTIA_MOTOR1S = 2.39e-4;		//!< [kgm^2] 1S軸 モータ側慣性
		static constexpr double INERTIA_MOTOR2L = 2.07e-4;		//!< [kgm^2] 2L軸 モータ側慣性
		static constexpr double INERTIA_MOTOR3U = 3.00e-5;		//!< [kgm^2] 3U軸 モータ側慣性
		static constexpr double VISCOSITY_ALL1S = 1.26e-2;		//!< [Nm/(rad/s)] 1S軸 全粘性 Da
		static constexpr double VISCOSITY_ALL2L = 9.30e-3;		//!< [Nm/(rad/s)] 2L軸 全粘性 Da
		static constexpr double VISCOSITY_ALL3U = 2.30e-3;		//!< [Nm/(rad/s)] 3U軸 全粘性 Da
		static constexpr double VISCOSITY_LOAD1S = 1.05;		//!< [Nm/(rad/s)] 1S軸 負荷側粘性
		static constexpr double VISCOSITY_LOAD2L = 12.0;		//!< [Nm/(rad/s)] 2L軸 負荷側粘性
		static constexpr double VISCOSITY_LOAD3U = 0.287;		//!< [Nm/(rad/s)] 3U軸 負荷側粘性
		static constexpr double VISCOSITY_MOTOR1S = 1.25e-2;	//!< [Nm/(rad/s)] 1S軸 モータ側粘性
		static constexpr double VISCOSITY_MOTOR2L = 8.10e-3;	//!< [Nm/(rad/s)] 2L軸 モータ側粘性
		static constexpr double VISCOSITY_MOTOR3U = 2.27e-3;	//!< [Nm/(rad/s)] 3U軸 モータ側粘性
		static constexpr double TORSION_VISCOSITY1S = 0.705;	//!< [Nm/(rad/s)] 1S軸 ねじれ粘性
		static constexpr double TORSION_VISCOSITY2L = 18.5;		//!< [Nm/(rad/s)] 2L軸 ねじれ粘性
		static constexpr double TORSION_VISCOSITY3U = 1.42;		//!< [Nm/(rad/s)] 3U軸 ねじれ粘性
		static constexpr double TORSION_STIFFNESS1S = 5.62e3;	//!< [Nm/rad] 1S軸 ねじれ剛性
		static constexpr double TORSION_STIFFNESS2L = 1.12e4;	//!< [Nm/rad] 2L軸 ねじれ剛性
		static constexpr double TORSION_STIFFNESS3U = 6.03e3;	//!< [Nm/rad] 3U軸 ねじれ剛性
		static constexpr double FRICTION_CLMB1S = 1.62e-1;		//!< [Nm] 1S軸 クーロン摩擦力 Fclmb
		static constexpr double FRICTION_CLMB2L = 1.54e-1;		//!< [Nm] 2L軸 クーロン摩擦力 Fclmb
		static constexpr double FRICTION_CLMB3U = 2.47e-2;		//!< [Nm] 3U軸 クーロン摩擦力 Fclmb
		static constexpr double HOMEPOS_X = 0.084166920164635;	//!< [m] ホームポジション作業空間位置 X軸 先端G(グラインド点)
		static constexpr double HOMEPOS_Y = 0;					//!< [m] ホームポジション作業空間位置 Y軸 先端G(グラインド点)
		static constexpr double HOMEPOS_Z = 0.365902851452706;	//!< [m] ホームポジション作業空間位置 Z軸 先端G(グラインド点)
		static constexpr double STBYPOS_X = 0.3000;				//!< [m] 待機ポジション作業空間位置 X軸 先端G(グラインド点)
		static constexpr double STBYPOS_Y = 0;					//!< [m] 待機ポジション作業空間位置 Y軸 先端G(グラインド点)
		static constexpr double STBYPOS_Z = 0.5250;				//!< [m] 待機ポジション作業空間位置 Z軸 先端G(グラインド点)
		static constexpr double FSENOFST_X = -24.1908;			//!< [N] ホームポジションでの力覚センサオフセット量 X軸
		static constexpr double FSENOFST_Y =  -0.4293;			//!< [N] ホームポジションでの力覚センサオフセット量 Y軸
		static constexpr double FSENOFST_Z =  24.1382;			//!< [N] ホームポジションでの力覚センサオフセット量 Z軸
		static constexpr double FSENGRAV_Z  = 33.9845;			//!< [N] 力覚センサから先端に掛かるツールの重力(作業空間Z軸)
		
		// ZFA01重力項の定義
		static constexpr double Ag3 = -17.5;	//!< [Nm]  重力項振幅定数 3U軸
		static constexpr double Pg3 = 0.785;	//!< [rad] 重力項位相定数 3U軸
		
		// 作業空間上の障害物の定義
		static constexpr double OBSTACLE1_Z =  50e-3;	//!< [m] 障害物1のZ方向位置の定義（床面）
		static constexpr double OBSTACLE2_Y = 200e-3;	//!< [m] 障害物2のY方向位置の定義（壁面）
		
		//! @brief トルク定数 [Nm/A]
		static constexpr Matrix<1,N_AX> Kt = {TORQUE_CONST1S, TORQUE_CONST2L, TORQUE_CONST3U};
		
		//! @brief 瞬時最大電流 [A]
		static constexpr Matrix<1,N_AX> iq_max = {CURRENT_MAX1S*0.5, CURRENT_MAX2L*0.5, CURRENT_MAX3U*0.5};
		
		//! @brief 緊急停止最大電流 [A]
		static constexpr Matrix<1,N_AX> iq_emer = {CURRENT_MAX1S*1.1, CURRENT_MAX2L*1.1, CURRENT_MAX3U*1.1};
		
		//! @brief 減速比 [-]
		static constexpr Matrix<1,N_AX> Rg = {GEAR_RATIO1S, GEAR_RATIO2L, GEAR_RATIO3U};
		
		//! @brief 全慣性 [kgm^2]
		static constexpr Matrix<1,N_AX> Ja = {INERTIA_ALL1S, INERTIA_ALL2L, INERTIA_ALL3U};
		
		//! @brief 全粘性 [Nm/(rad/s)]
		static constexpr Matrix<1,N_AX> Da = {VISCOSITY_ALL1S, VISCOSITY_ALL2L, VISCOSITY_ALL3U};
		
		//! @brief 負荷側慣性 [kgm^2]
		static constexpr Matrix<1,N_AX> Jl = {INERTIA_LOAD1S, INERTIA_LOAD2L, INERTIA_LOAD3U};
		
		//! @brief 負荷側粘性 [Nm/(rad/s)]
		static constexpr Matrix<1,N_AX> Dl = {VISCOSITY_LOAD1S, VISCOSITY_LOAD2L, VISCOSITY_LOAD3U};
		
		//! @brief ねじれ粘性 [Nm/(rad/s)]
		static constexpr Matrix<1,N_AX> Ds = {TORSION_VISCOSITY1S, TORSION_VISCOSITY2L, TORSION_VISCOSITY3U};
		
		//! @brief ねじれ剛性 [Nm/rad]
		static constexpr Matrix<1,N_AX> Ks = {TORSION_STIFFNESS1S, TORSION_STIFFNESS2L, TORSION_STIFFNESS3U};
		
		//! @brief モータ側慣性 [kgm^2]
		static constexpr Matrix<1,N_AX> Jm = {INERTIA_MOTOR1S, INERTIA_MOTOR2L, INERTIA_MOTOR3U};
		
		//! @brief モータ側粘性 [Nm/(rad/s)]
		static constexpr Matrix<1,N_AX> Dm = {VISCOSITY_MOTOR1S, VISCOSITY_MOTOR2L, VISCOSITY_MOTOR3U};
		
		//! @brief クーロン摩擦力 [Nm]
		static constexpr Matrix<1,N_AX> Fclmb = {FRICTION_CLMB1S, FRICTION_CLMB2L, FRICTION_CLMB3U};
		
		//! @brief 2慣性系パラメータ構造体 1S軸
		static constexpr struct TwoInertiaParams Params1S = {Jl[1], Dl[1], Ds[1], Ks[1], Jm[1], Dm[1], Rg[1], Kt[1]};
		
		//! @brief 2慣性系パラメータ構造体 2L軸
		static constexpr struct TwoInertiaParams Params2L = {Jl[2], Dl[2], Ds[2], Ks[2], Jm[2], Dm[2], Rg[2], Kt[2]};
		
		//! @brief 2慣性系パラメータ構造体 3U軸
		static constexpr struct TwoInertiaParams Params3U = {Jl[3], Dl[3], Ds[3], Ks[3], Jm[3], Dm[3], Rg[3], Kt[3]};
		
		//! @brief 2慣性系パラメータ構造体 全軸
		static constexpr std::array<struct TwoInertiaParams, N_AX> ParamsAll = {Params1S, Params2L, Params3U};
		
		//! @brief ホームポジション作業空間位置ベクトル [m] 先端G(グラインド点)
		static constexpr Matrix<1,3> pGhome = {HOMEPOS_X, HOMEPOS_Y, HOMEPOS_Z};
		
		//! @brief 待機ポジション作業空間位置ベクトル [m] 先端G(グラインド点)
		static constexpr Matrix<1,3> pGstby = {STBYPOS_X, STBYPOS_Y, STBYPOS_Z};
		
		//! @brief 待機ポジションでの力覚センサオフセット量 [N]
		static constexpr Matrix<1,3> fofst = {FSENOFST_X, FSENOFST_Y, FSENOFST_Z};
		
		//! @brief 力覚センサから先端に掛かるツールの重力 [N]
		static constexpr Matrix<1,3> fgrav = {0, 0, FSENGRAV_Z};
		
		ZFA01();			//!< コンストラクタ
		ZFA01(ZFA01&& r);	//!< ムーブコンストラクタ
		~ZFA01();			//!< デストラクタ
		
		void CurrentLimiter(Matrix<1,N_AX>& CurrentRef);		//!< 電流指令リミッタ
		void ConvJointAngleToMDHparam(const Matrix<1,N_AX>& JointAngle, Matrix<1,N_AX>& JointAngleMDH) const;	//!< 関節角から修正DHパラメータ定義上の関節角に換算する関数
		void ConvImplementVecToMDHparam(const Matrix<1,N_AX>& ImplementVec, Matrix<1,N_AX>& MDHVec) const;		//!< 実装上のベクトルから修正DHパラメータ定義上のベクトルに換算する関数
		void ConvMDHparamToImplementVec(const Matrix<1,N_AX>& MDHVec, Matrix<1,N_AX>& ImplementVec) const;		//!< 修正DHパラメータ定義上のベクトルから実装上のベクトルに換算する関数
		void GetWorkspacePosition(const Matrix<1,N_AX>& JointAngle, std::array<Matrix<1,N_AX>, 4>& WorkPos);	//!< 作業空間位置を順運動学により計算する関数(全軸分を計算する版)
		void GetWorkspacePosition(const Matrix<1,N_AX>& JointAngle, Matrix<1,N_AX>& p1, Matrix<1,N_AX>& p2, Matrix<1,N_AX>& p3, Matrix<1,N_AX>& pG);	//!< 作業空間位置を順運動学により計算する関数(全軸分を計算する版, 個別に返す版)
		void GetWorkspaceForce(const Matrix<1,3>& SensForce, Matrix<1,3>& fG);	//! @brief 作業空間並進推力を順運動学により計算する関数(個別に返す版) 注意！： GetWorkspacePosition()が予め実行済みであること。
		void UpdateJacobi(const Matrix<1,N_AX>& JointAngle);			//!< Jacobi行列を更新する関数
		Matrix<1,N_AX> Jaco(const Matrix<1,N_AX>& JointSpace) const;	//!< Jacobi行列とベクトルの積を計算してベクトルを返す関数
		Matrix<1,N_AX> Jacoinv(const Matrix<1,N_AX>& WorkSpace) const;	//!< 逆Jacobi行列とベクトルの積を計算してベクトルを返す関数
		Matrix<1,N_AX> JacoT(const Matrix<1,N_AX>& WorkSpace) const;	//!< 転置Jacobi行列とベクトルの積を計算してベクトルを返す関数
		Matrix<1,N_AX> JacoTinv(const Matrix<1,N_AX>& JointSpace) const;//!< 転置の逆Jacobi行列とベクトルの積を計算してベクトルを返す関数
		void CheckAngleSafety(const Matrix<1,N_AX>& JointAngle) const;	//!< 関節可動範囲の安全確認(範囲外になると緊急停止)
		void CheckRangeSafety(const Matrix<1,N_AX>& WorkspacePos) const;//!< 作業可動範囲の安全確認(範囲外になると緊急停止)
		void CheckSafety(const Matrix<1,N_AX>& JointAngle, const Matrix<1,N_AX>& WorkspacePos) const;	//!< ロボットの安全確認（範囲外になると緊急停止）
		Matrix<1,N_AX> GetGravityCurrent(const Matrix<1,N_AX>& JointAngle) const;		//!< 重力項分の補償電流を取得する関数
		Matrix<1,N_AX> GetGravityTorque(const Matrix<1,N_AX>& JointAngle) const;		//!< 重力項分の補償トルクを取得する関数
		
	private:
		ZFA01(const ZFA01&) = delete;					//!< コピーコンストラクタ使用禁止
		const ZFA01& operator=(const ZFA01&) = delete;	//!< 代入演算子使用禁止
		
		// 修正Denavit-Hartenbergパラメータの定義(ZFA01のCADデータから導出)
		static constexpr double Rdisk = 103e-3;	//!< [m] ディスクグラインダの直径
		static constexpr double la = 406e-3;	//!< [m] リンクオフセット
		static constexpr double lb = 300e-3;	//!< [m] リンクオフセット
		static constexpr double lc =  50e-3;	//!< [m] リンクオフセット
		static constexpr double ld = 300e-3;	//!< [m] リンクオフセット
		static constexpr double le = 179.47e-3 + Rdisk/2.0; //!< [m] リンクオフセット
		
		// 各軸間の同次変換行列
		Matrix<4,4> T01;	//!< ベース0から1S軸までの同次変換行列
		Matrix<4,4> T12;	//!< 1S軸から2L軸までの同次変換行列
		Matrix<4,4> T23;	//!< 2L軸から3U軸までの同次変換行列
		Matrix<4,4> T3g;	//!< 3U軸から先端G(グラインド点)までの同次変換行列
		Matrix<4,4> T02;	//!< ベース0から2L軸までの同次変換行列
		Matrix<4,4> T03;	//!< ベース0から3U軸までの同次変換行列
		Matrix<4,4> T0g;	//!< ベース0から先端G(グラインド点)までの同次変換行列
		
		// Jacobi行列の群れ
		Matrix<N_AX,N_AX> J;	//!< Jacobi行列
		Matrix<N_AX,N_AX> Jinv;	//!< 逆Jacobi行列
		Matrix<N_AX,N_AX> Jt;	//!< 転置Jacobi行列
		Matrix<N_AX,N_AX> Jtinv;//!< 転置の逆Jacobi行列
};
}

#endif

