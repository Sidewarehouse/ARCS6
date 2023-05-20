//! @file ResoLoopCutDistObsrv.hh
//! @brief 共振ループ切断外乱オブザーバ
//!
//! モータ側トルク指令とねじれトルクセンサ値から，モータ側外乱と負荷側速度外乱を推定するオブザーバで，
//! 2慣性共振系の共振ループを切断するための特殊な外乱オブザーバのクラス
//!
//! @date 2022/03/30
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RESOLOOPCUTDISTOBSRV
#define RESOLOOPCUTDISTOBSRV

#include <cassert>
#include "Matrix.hh"
#include "TransferFunction.hh"
#include "TwoInertiaParamDef.hh"

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

namespace ARCS {	// ARCS名前空間
//! @brief 共振ループ切断外乱オブザーバ
//! @tparam 
//template <>
class ResoLoopCutDistObsrv {
	public:
		//! @brief 空コンストラクタ
		ResoLoopCutDistObsrv(void)
			: Ds(), Ks(), Jm(), Dm(), Rg(), gd(), Ts(),
			  Lnum(), Lden(), Rnum(), Rden(), Qnum(), Qden(), L(), R(), QLinv()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ
		//! @param[in]	Params		2慣性共振系パラメータ構造体
		//! @param[in]	Bandwidth	推定帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		ResoLoopCutDistObsrv(const struct TwoInertiaParams& Params, const double Bandwidth, const double SmplTime)
			: Ds(), Ks(), Jm(), Dm(), Rg(), gd(), Ts(),
			  Lnum(), Lden(), Rnum(), Rden(), Qnum(), Qden(), L(), R(), QLinv()
		{
			SetParameters(Params, Bandwidth, SmplTime);	// オブザーバのパラメータを設定
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		ResoLoopCutDistObsrv(ResoLoopCutDistObsrv&& r)
			: Ds(r.Ds), Ks(r.Ks), Jm(r.Jm), Dm(r.Dm), Rg(r.Rg), gd(r.gd), Ts(r.Ts),
			  Lnum(r.Lnum), Lden(r.Lden), Rnum(r.Rnum), Rden(r.Rden), Qnum(r.Qnum), Qden(r.Qden), L(), R(), QLinv()
		{
			
		}
		
		//! @brief デストラクタ
		~ResoLoopCutDistObsrv(){
			PassedLog();
		}
		
		//! @brief オブザーバのパラメータを設定する関数
		//! @param[in]	Params		2慣性共振系パラメータ構造体
		//! @param[in]	Bandwidth	推定帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		void SetParameters(const struct TwoInertiaParams& Params, const double Bandwidth, const double SmplTime){
			// パラメータの格納
			Ds = Params.Ds;
			Ks = Params.Ks;
			Jm = Params.Jm;
			Dm = Params.Dm;
			Rg = Params.Rg;
			gd = Bandwidth;
			Ts = SmplTime;
			
			// オブザーバ左側ループの伝達関数
			Lnum.Set( 1.0/Rg );
			Lden.Set( Jm, Dm );
			L.SetCoefficients(Lnum, Lden, Ts);
			
			// オブザーバ右側ループの伝達関数
			Rnum.Set( 1, 0 );
			Rden.Set( Ds, Ks );
			R.SetCoefficients(Rnum, Rden, Ts);
			
			// オブザーバフィルタと左側ループの逆系の乗算結果の伝達関数
			Qnum.Set( Rg*gd*Jm, Rg*gd*Dm );
			Qden.Set( 1, gd );
			QLinv.SetCoefficients(Qnum, Qden, Ts);
		}
		
		//! @brief オブザーバの制御帯域を設定する関数
		//! @param[in]	Bandwidth	推定帯域 [rad/s]
		void SetBandwidth(const double Bandwidth){
			// パラメータの格納
			gd = Bandwidth;
			
			// オブザーバフィルタと左側ループの逆系の乗算結果の伝達関数
			Qnum.Set( Rg*gd*Jm, Rg*gd*Dm );
			Qden.Set( 1, gd );
			QLinv.SetCoefficients(Qnum, Qden, Ts);
		}
		
		//! @brief モータ側トルク補償値を取得する関数
		//! @param[in]	TorqueRef		モータ側トルク指令 [Nm]
		//! @param[in]	TorsionTorque	ねじれトルク [Nm]
		//! @return	モータ側トルク補償値 [Nm]
		double GetCompTorque(const double TorqueRef, const double TorsionTorque){
			return QLinv.GetResponse( L.GetResponse(TorqueRef) - R.GetResponse(TorsionTorque) );	// [Nm] モータ側トルク補償値の計算
		}
		
		//! @brief RLC-TTC用のP-D制御ゲインを計算する関数(仮設)
		//! @param[in]	wt	ねじれトルク制御帯域 [rad/s]
		//! @param[in]	zt	ねじれトルク制御制動係数 [-]
		//! @param[out]	Kpt	P-Dねじれトルク制御器 比例ゲイン
		//! @param[out] Kdt P-Dねじれトルク制御器 微分ゲイン
		void GetPDgainForRLCTTC(const double wt, const double zt, double& Kpt, double& Kdt){
			Kpt = -( (Jm*Ks - Dm*Ds)*Rg*wt*wt )/( 2.0*Ds*Ks*wt*zt - Ds*Ds*wt*wt - Ks*Ks );						// [-] P-Dねじれトルク制御 比例ゲイン
			Kdt = -( 2.0*Jm*Ks*Rg*wt*zt - Ds*Jm*Rg*wt*wt - Dm*Ks*Rg )/( 2.0*Ds*Ks*wt*zt - Ds*Ds*wt*wt - Ks*Ks );// [-] P-Dねじれトルク制御 微分ゲイン
		}
		
	private:
		ResoLoopCutDistObsrv(const ResoLoopCutDistObsrv&) = delete;					//!< コピーコンストラクタ使用禁止
		const ResoLoopCutDistObsrv& operator=(const ResoLoopCutDistObsrv&) = delete;//!< 代入演算子使用禁止
		double Ds;	//!< [Nm/(rad/s)] ねじれ粘性
		double Ks;	//!< [Nm/rad] ねじれ剛性
		double Jm;	//!< [kgm^2] モータ側慣性
		double Dm;	//!< [Nm/(rad/s)] モータ側粘性
		double Rg;	//!< [-] 減速比
		double gd;	//!< [rad/s] 推定帯域
		double Ts;	//!< [s] サンプリング周期
		Matrix<1,1> Lnum;	//!< オブザーバ左側ループの伝達関数の分子係数ベクトル
		Matrix<1,2> Lden;	//!< オブザーバ左側ループの伝達関数の分母係数ベクトル
		Matrix<1,2> Rnum;	//!< オブザーバ右側ループの伝達関数の分子係数ベクトル
		Matrix<1,2> Rden;	//!< オブザーバ右側ループの伝達関数の分母係数ベクトル
		Matrix<1,2> Qnum;	//!< オブザーバフィルタと左側ループの逆系の乗算結果の伝達関数の分子係数ベクトル
		Matrix<1,2> Qden;	//!< オブザーバフィルタと左側ループの逆系の乗算結果の伝達関数の分母係数ベクトル
		TransferFunction<0,1> L;	//!< 外乱オブザーバの左側の伝達関数
		TransferFunction<1,1> R;	//!< 外乱オブザーバの右側の伝達関数
		TransferFunction<1,1> QLinv;//!< 外乱オブザーバのLPFと左側の逆系の乗算結果の伝達関数
};
}

#endif

