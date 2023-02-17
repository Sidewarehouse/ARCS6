//! @file ZFA01.cc
//! @brief Zero Force Arm 01 - ZFA01クラス
//!
//! ZFA01ロボット制御のためのクラス
//!
//! @date 2022/03/18
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include "ZFA01.hh"
#include "Limiter.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

using namespace ARCS;

//! @brief コンストラクタ
ZFA01::ZFA01()
	: T01(), T12(), T23(), T3g(), T02(), T03(), T0g(),
	  J(), Jinv(), Jt(), Jtinv()
{
	PassedLog();
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
ZFA01::ZFA01(ZFA01&& r)
	: T01(r.T01), T12(r.T12), T23(r.T23), T3g(r.T3g), T02(r.T02), T03(r.T03), T0g(r.T0g),
	  J(r.J), Jinv(r.Jinv), Jt(r.Jt), Jtinv(r.Jtinv)
{
	
}

//! @brief デストラクタ
ZFA01::~ZFA01(){
	PassedLog();
}

//! @brief 電流指令リミッタ
//! @param[in,out]	CurrentRef	電流指令ベクトル [A]
void ZFA01::CurrentLimiter(Matrix<1,N_AX>& CurrentRef){
	arcs_assert(CurrentRef[1] < std::abs(ZFA01::iq_emer[1]));	// [A] 緊急停止最大電流 1S軸
	arcs_assert(CurrentRef[2] < std::abs(ZFA01::iq_emer[2]));	// [A] 緊急停止最大電流 2L軸
	arcs_assert(CurrentRef[3] < std::abs(ZFA01::iq_emer[3]));	// [A] 緊急停止最大電流 3U軸
	Limiter<ZFA01::N_AX>(CurrentRef, ZFA01::iq_max);			// [A] q軸電流リミッタ
}

//! @brief 関節角から修正DHパラメータ定義上の関節角に換算する関数
//! @param[in]	JointAngle		[rad] 関節角
//! @param[out]	JointAngleMDH	[rad] 修正DHパラメータ定義上の関節角
void ZFA01::ConvJointAngleToMDHparam(const Matrix<1,N_AX>& JointAngle, Matrix<1,N_AX>& JointAngleMDH) const{
	// ホームポジション(エンコーダ原点)から修正DHパラメータ定義上の原点への変換
	// ロボットに貼ってあるシールの回転方向定義から修正DHパラメータ定義の回転方向へも同時に変換
	JointAngleMDH[1] =  JointAngle[1];				// [rad] 1S軸はそのまま
	JointAngleMDH[2] =  JointAngle[2] + M_PI_2;		// [rad] 2L軸は正回転で，＋90度オフセット
	JointAngleMDH[3] = -JointAngle[3] - M_PI/4.0;	// [rad] 3U軸は逆回転で，ー45度オフセット
}

//! @brief 実装上のベクトルから修正DHパラメータ定義上のベクトルに換算する関数
//! @param[in]	ImplementVec	[*] 実装上のベクトル
//! @param[out]	MDHVec			[*] 修正DHパラメータ定義上のベクトル
void ZFA01::ConvImplementVecToMDHparam(const Matrix<1,N_AX>& ImplementVec, Matrix<1,N_AX>& MDHVec) const{
	// ロボットに貼ってあるシールの回転方向定義に従って符号を変換
	MDHVec[1] =  ImplementVec[1];	// [*] 1S軸はそのまま
	MDHVec[2] =  ImplementVec[2];	// [*] 2L軸もそのまま
	MDHVec[3] = -ImplementVec[3];	// [*] 3U軸は逆回転
}

//! @brief 修正DHパラメータ定義上のベクトルから実装上のベクトルに換算する関数
//! @param[in]	MDHVec			[*] 修正DHパラメータ定義上のベクトル
//! @param[out]	ImplementVec	[*] 実装上のベクトル
void ZFA01::ConvMDHparamToImplementVec(const Matrix<1,N_AX>& MDHVec, Matrix<1,N_AX>& ImplementVec) const{
	// ロボットに貼ってあるシールの回転方向定義に従って符号を変換
	ImplementVec[1] =  MDHVec[1];	// [*] 1S軸はそのまま
	ImplementVec[2] =  MDHVec[2];	// [*] 2L軸もそのまま
	ImplementVec[3] = -MDHVec[3];	// [*] 3U軸は逆回転
}

//! @brief 作業空間位置を順運動学により計算する関数(全軸分を計算する版, 配列で返す版)
//! @param[in]	JointAngle	 [rad] 3軸関節角度ベクトル [th1,th2,th3]^T
//! @param[out]	WorkPosition [m] 1S軸,2L軸,3U軸,先端G(グラインド点)の作業空間位置の配列
void ZFA01::GetWorkspacePosition(const Matrix<1,N_AX>& JointAngle, std::array<Matrix<1,N_AX>, 4>& WorkPosition){
	// 個別に返す版で計算して配列に格納
	GetWorkspacePosition(JointAngle, WorkPosition.at(0), WorkPosition.at(1), WorkPosition.at(2), WorkPosition.at(3));
}

//! @brief 作業空間位置を順運動学により計算する関数(全軸分を計算する版, 個別に返す版)
//! @param[in]	JointAngle	 [rad] 3軸関節角度ベクトル [th1,th2,th3]^T
//! @param[out]	p1	[m] 1S軸の作業空間位置
//! @param[out]	p2	[m] 2L軸の作業空間位置
//! @param[out]	p3	[m] 3U軸の作業空間位置
//! @param[out]	pG	[m] 先端G(グラインド点)の作業空間位置
void ZFA01::GetWorkspacePosition(const Matrix<1,N_AX>& JointAngle, Matrix<1,N_AX>& p1, Matrix<1,N_AX>& p2, Matrix<1,N_AX>& p3, Matrix<1,N_AX>& pG){
	Matrix<1,N_AX> theta;						// [rad] 修正DHパラメータ定義上の関節角
	ConvJointAngleToMDHparam(JointAngle, theta);// 修正DHパラメータ定義上の関節角に変換
	
	// 正弦と余弦の予備計算
	const double S1 = sin(theta[1]);
	const double C1 = cos(theta[1]);
	const double S2 = sin(theta[2]);
	const double C2 = cos(theta[2]);
	const double S3 = sin(theta[3]);
	const double C3 = cos(theta[3]);
	
	// 各軸間の同次変換行列
	T01.Set(
		C1, -S1,  0,  0,
		S1,  C1,  0,  0,
		 0,   0,  1, la,
		 0,   0,  0,  1
	);
	T12.Set(
		C2, -S2,  0,  0,
		 0,   0, -1,  0,
		S2,  C2,  0,  0,
		 0,   0,  0,  1
	);
	T23.Set(
		C3, -S3,  0, lb,
		S3,  C3,  0,  0,
		 0,   0,  1,  0,
		 0,   0,  0,  1
	);
	T3g.Set(
		 1,  0,  0, lc-le,
		 0,  0, -1,   -ld,
		 0,  1,  0,     0,
		 0,  0,  0,     1
	);
	
	// ベース0から各軸への同次変換行列
	T02 = T01*T12;
	T03 = T02*T23;
	T0g = T03*T3g;
	
	// 各々の同次変換行列から作業空間位置を抽出
	for(size_t i = 1; i <= 3; ++i){
		p1[i] = T01.GetElem(i,4);	// [m] 1S軸の作業空間位置を抽出
		p2[i] = T02.GetElem(i,4);	// [m] 2L軸の作業空間位置を抽出
		p3[i] = T03.GetElem(i,4);	// [m] 3U軸の作業空間位置を抽出
		pG[i] = T0g.GetElem(i,4);	// [m] 先端G(グラインド点)の作業空間位置を抽出
	}
}

//! @brief 作業空間並進推力を順運動学により計算する関数(個別に返す版)
//!        注意！： GetWorkspacePosition()が予め実行済みであること。
//! @param[in]	SensForce	[N] 力覚センサの検出推力ベクトル [fsen_x, fsen_y, fsen_z]^T
//! @param[out]	fG			[N] 先端G(グラインド点)の作業空間推力ベクトル [fG_x, fG_y, fG_z]^T
void ZFA01::GetWorkspaceForce(const Matrix<1,3>& SensForce, Matrix<1,3>& fG){
	Matrix<3,3> R0g;						// 先端Gからベース0への回転行列
	getsubmatrix(T0g, 1, 1, R0g);			// 同次変換行列から回転行列を抽出
	fG = R0g*(SensForce - fofst) - fgrav;	// [N] 先端Gの作業空間推力に変換
}

//! @brief Jacobi行列を更新する関数
//! @param[in]	JointAngle	3軸関節角度ベクトル [th1,th2,th3]^T [rad]
void ZFA01::UpdateJacobi(const Matrix<1,N_AX>& JointAngle){
	Matrix<1,N_AX> theta;						// [rad] 修正DHパラメータ定義上の関節角
	ConvJointAngleToMDHparam(JointAngle, theta);// 修正DHパラメータ定義上の関節角に変換
	
	// 正弦と余弦の計算
	const double S1 = sin(theta[1]);
	const double C1 = cos(theta[1]);
	const double S2 = sin(theta[2]);
	const double C2 = cos(theta[2]);
	const double S3 = sin(theta[3]);
	const double C3 = cos(theta[3]);
	const double S23 = sin(theta[2] + theta[3]);
	const double C23 = cos(theta[2] + theta[3]);
	
	// Jacobi行列 (ZFA01G10VE_03_JacobiMat.macの計算結果)
	J.Set(
		 ((le-lc)*C23-lb*C2)*S1-ld*S1*S23, ((le-lc)*C1*C2-ld*C1*S2)*S3+((le-lc)*C1*C3-lb*C1)*S2+ld*C1*C2*C3, ((le-lc)*C1*C2-ld*C1*S2)*S3+(le-lc)*C1*C3*S2+ld*C1*C2*C3,
		ld*C1*S23+(lc-le)*C1*C23+lb*C1*C2,    ((le-lc)*C2*S1-ld*S1*S2)*S3+((le-lc)*C3-lb)*S1*S2+ld*C2*C3*S1, ((le-lc)*C2*S1-ld*S1*S2)*S3+(le-lc)*C3*S1*S2+ld*C2*C3*S1,
		                                0,               ((le-lc)*S2+ld*C2)*S3+ld*C3*S2+(lc-le)*C2*C3+lb*C2,             ((le-lc)*S2+ld*C2)*S3+ld*C3*S2+(lc-le)*C2*C3
	);
	Jt = tp(J);				// 転置Jacobi行列
	
	// 特異点に入って無ければ逆Jacobi行列を更新
	if(EPS_SINGULAR < std::abs(det(J))){
		Jinv = inv(J);		// 逆Jacobi行列
		Jtinv = inv(tp(J));	// 転置の逆Jacobi行列
	}
}

//! @brief Jacobi行列とベクトルの積を計算してベクトルを返す関数
//! @param[in]	JointSpace	関節空間ベクトル
//! @return	作業空間ベクトル
Matrix<1,ZFA01::N_AX> ZFA01::Jaco(const Matrix<1,N_AX>& JointSpace) const{
	Matrix<1,N_AX> v;
	ConvImplementVecToMDHparam(JointSpace, v);		// 符号を修正DHパラメータ定義上に合わせる
	return J*v;		// Jacobi行列を掛けて返す
}

//! @brief 逆Jacobi行列とベクトルの積を計算してベクトルを返す関数
//! @param[in]	WorkSpace	作業空間ベクトル
//! @return 関節空間ベクトル
Matrix<1,ZFA01::N_AX> ZFA01::Jacoinv(const Matrix<1,N_AX>& WorkSpace) const{
	Matrix<1,N_AX> v;
	ConvMDHparamToImplementVec(Jinv*WorkSpace, v);	// 符号をロボット回転方向定義に合わせる
	return v;
}

//! @brief 転置Jacobi行列とベクトルの積を計算してベクトルを返す関数
//! @param[in]	WorkSpace	作業空間ベクトル
//! @return 関節空間ベクトル
Matrix<1,ZFA01::N_AX> ZFA01::JacoT(const Matrix<1,N_AX>& WorkSpace) const{
	Matrix<1,N_AX> v;
	ConvMDHparamToImplementVec(Jt*WorkSpace, v);	// 符号をロボット回転方向定義に合わせる
	return v;
}

//! @brief 転置の逆Jacobi行列とベクトルの積を計算してベクトルを返す関数
//! @param[in]	JointSpace	関節空間ベクトル
//! @return	作業空間ベクトル
Matrix<1,ZFA01::N_AX> ZFA01::JacoTinv(const Matrix<1,N_AX>& JointSpace) const{
	Matrix<1,N_AX> v;
	ConvImplementVecToMDHparam(JointSpace, v);		// 符号を修正DHパラメータ定義上に合わせる
	return Jtinv*v;	// 転置の逆Jacobi行列を掛けて返す
}

//! @brief 関節可動範囲の安全確認（範囲外になると緊急停止）
//! @param[in]	JointAngle	[rad] 3軸関節角度ベクトル[th1,th2,th3]^T
void ZFA01::CheckAngleSafety(const Matrix<1,N_AX>& JointAngle) const{
	// 可動範囲監視
	arcs_assert(POSITION_LIMNEG1 < JointAngle[1] && JointAngle[1] < POSITION_LIMPOS1);	// 可動範囲チェック 1S軸
	arcs_assert(POSITION_LIMNEG2 < JointAngle[2] && JointAngle[2] < POSITION_LIMPOS2);	// 可動範囲チェック 2L軸
	arcs_assert(POSITION_LIMNEG3 < JointAngle[3] && JointAngle[3] < POSITION_LIMPOS3);	// 可動範囲チェック 3U軸
}

//! @brief 作業可動範囲の安全確認（範囲外になると緊急停止）
//! @param[in]	WorkspacePos	[m] 作業空間位置ベクトル[px,py,pz]^T
void ZFA01::CheckRangeSafety(const Matrix<1,N_AX>& WorkspacePos) const{
	// 可動範囲監視（不等式を満たしたとき正常）
	arcs_assert(OBSTACLE1_Z < WorkspacePos[3]);	// 可動範囲チェック 床面
	arcs_assert(WorkspacePos[2] < OBSTACLE2_Y);	// 可動範囲チェック 壁面
}

//! @brief ロボットの安全確認（範囲外になると緊急停止）
//! @param[in]	JointAngle		[rad] 3軸関節角度ベクトル[th1,th2,th3]^T
//! @param[in]	WorkspacePos	[m]   作業空間位置ベクトル[px,py,pz]^T 
void ZFA01::CheckSafety(const Matrix<1,N_AX>& JointAngle, const Matrix<1,N_AX>& WorkspacePos) const{
	// 安全監視
	CheckAngleSafety(JointAngle);	// 関節可動範囲のチェック
	CheckRangeSafety(WorkspacePos);	// 作業可動範囲のチェック
}

//! @brief 重力項分の補償電流を取得する関数
//! @param[in]	JointAngle	[rad] 3軸関節角度ベクトル[th1,th2,th3]^T
//! @return	[A] 重力補償電流
Matrix<1,ZFA01::N_AX> ZFA01::GetGravityCurrent(const Matrix<1,N_AX>& JointAngle) const{
	Matrix<1,N_AX> igrav;
	igrav[1] = 0;	// [A] 重力補償電流 1S軸
	igrav[2] = 0;	// [A] 重力補償電流 2L軸
	igrav[3] = Ag3*cos(JointAngle[3] + Pg3)/(ZFA01::Kt[3]*ZFA01::Rg[3]);	// [A] 重力補償電流 3U軸
	return igrav;
}

//! @brief 重力項分の補償トルクを取得する関数
//! @param[in]	JointAngle	[rad] 3軸関節角度ベクトル[th1,th2,th3] ^T
//! @return	[Nm] 重力補償トルク
Matrix<1,ZFA01::N_AX> ZFA01::GetGravityTorque(const Matrix<1,N_AX>& JointAngle) const{
	Matrix<1,N_AX> tgrav;
	tgrav[1] = 0;	// [A] 重力補償電流 1S軸
	tgrav[2] = 0;	// [A] 重力補償電流 2L軸
	tgrav[3] = Ag3*cos(JointAngle[3] + Pg3);	// [Nm] 重力補償トルク 3U軸
	return tgrav;
}
