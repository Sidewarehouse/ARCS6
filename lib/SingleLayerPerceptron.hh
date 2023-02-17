//! @file SingleLayerPerceptron.hh
//! @brief 単層パーセプトロンクラス
//!
//! 1層のパーセプトロンレイヤ
//!
//! @date 2021/12/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SINGLELAYERPERCEPTRON
#define SINGLELAYERPERCEPTRON

#include <cassert>
#include <array>
#include <cmath>
#include <string>
#include "Matrix.hh"
#include "NeuralNetParamDef.hh"
#include "ActivationFunctions.hh"
#include "RandomGenerator.hh"
#include "Statistics.hh"
#include "CsvManipulator.hh"

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
//! @brief 単層パーセプトロンクラス
//! @tparam N	入力の数
//! @tparam P	パーセプトロンの数
//! @tparam M	ミニバッチ数
//! @tparam	AF	活性化関数の種類
//! @tparam	IT	初期化のタイプ
//! @tparam	GD	勾配降下法のタイプ
//! @tparam	DD	ドロップアウト動作フラグ
template <
	size_t N,
	size_t P,
	size_t M,
	ActvFunc AF,
	NnInitTypes IT = NnInitTypes::XAVIER,
	NnDescentTypes GD = NnDescentTypes::MOMENTUM,
	NnDropout DD = NnDropout::DISABLE
>
class SingleLayerPerceptron {
	public:
		//! @brief コンストラクタ
		SingleLayerPerceptron()
			: u(), U(), W(), b(), l(Matrix<1,M>::ones()), lT(Matrix<M,1>::ones()),
				delta(), Delta(), fpu(), fpU(), dW(), db(), DW(), Db(),
				HW(), Hb(), GW(), Gb(),
				HWhat(), GWhat(), Hbhat(), Gbhat(), SGDcount(0),
				DropRand(0, 1), DropMask(),
				xbar(0), sigma(1), eps(0.01), alph(0.9), bet(0.999), NearZero(1e-8), DropRate(0.5)
		{
			
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		SingleLayerPerceptron(SingleLayerPerceptron&& r)
			: u(r.u), U(r.U), W(r.W), b(r.b), l(r.l), lT(r.lT),
				delta(r.delta), Delta(r.Delta), fpu(r.fpu), fpU(r.fpU), dW(r.dW), db(r.db),
				HW(r.HW), Hb(r.Hb), GW(r.GW), Gb(r.Gb),
				HWhat(r.HWhat), GWhat(r.GWhat), Hbhat(r.Hbhat), Gbhat(r.Gbhat), SGDcount(r.SGDcount),
				DropRand(r.DropRand), DropMask(r.DropMask),
				xbar(r.xbar), sigma(r.sigma), eps(r.eps), alph(r.alph), bet(r.bet), NearZero(r.NearZero), DropRate(r.DropRate)
		{
			
		}
		
		//! @brief デストラクタ
		~SingleLayerPerceptron(){
			
		}
		
		//! @brief 重み行列の正規分布乱数による初期化
		//! @param	sigma	正規分布乱数の標準偏差
		void InitWeightUsingGaussianRandom(const double sigma){
			RandomGenerator Rand(0, sigma);	// メルセンヌ・ツイスタの生成
			Rand.GetGaussianRandomMatrix(W);// 平均0，標準偏差σのガウシアン乱数行列の取得
		}
		
		//! @brief 重み行列の初期化
		//! @param[in]	Nprev	前層のユニット数
		void InitWeight(size_t Nprev){
			if constexpr(IT == NnInitTypes::XAVIER){
				// Xavierの初期化
				InitWeightUsingGaussianRandom( 1.0/sqrt((double)Nprev) );
			}
			if constexpr(IT == NnInitTypes::HE){
				// Heの初期化
				InitWeightUsingGaussianRandom( sqrt(2.0/(double)Nprev) );
			}
		}
		
		//! @brief 確率的勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	epsilon	更新ゲイン(学習率)
		void SetGainOfSGD(double epsilon){
			eps = epsilon;
		}
		
		//! @brief モーメンタム確率的勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	epsilon	更新ゲイン(学習率)
		//! @param[in]	alpha	更新ゲイン(運動量項の学習率)
		void SetGainOfMomentumSGD(double epsilon, double alpha){
			eps = epsilon;
			alph = alpha;
		}
		
		//! @brief AdaGrad勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	epsilon	更新ゲイン(学習率)
		//! @param[in]	zero	ゼロ割回避用の係数
		void SetGainOfAdaGrad(double epsilon, double zero){
			eps = epsilon;
			NearZero = zero;
		}
		
		//! @brief RMSprop勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	epsilon	更新ゲイン(学習率)
		//! @param[in]	alpha	更新ゲイン(学習率)
		//! @param[in]	zero	ゼロ割回避用の係数
		void SetGainOfRMSprop(double epsilon, double alpha, double zero){
			eps = epsilon;
			alph = alpha;
			NearZero = zero;
		}

		//! @brief AdaDelta勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	alpha	更新ゲイン(学習率)
		//! @param[in]	zero	ゼロ割回避用の係数
		void SetGainOfAdaDelta(double alpha, double zero){
			alph = alpha;
			NearZero = zero;
		}
		
		//! @brief Adam勾配降下法の更新ゲイン(学習率)の設定
		//! @param[in]	epsilon	更新ゲイン(学習率)
		//! @param[in]	alpha	更新ゲイン(学習率)
		//! @param[in]	beta	更新ゲイン(学習率)
		//! @param[in]	zero	ゼロ割回避用の係数
		void SetGainOfAdam(double epsilon, double alpha, double beta, double zero){
			eps = epsilon;
			alph = alpha;
			bet = beta;
			NearZero = zero;
			SGDcount = 0;
		}
		
		//! @brief ドロップアウト率の設定
		//! @param[in]	DropoutRate	ドロップアウト率(1のときドロップアウトしない，0.5のとき半分ドロップアウト，0のとき全部ドロップアウト)
		void SetDropoutRate(double DropoutRate){
			DropRate = DropoutRate;
		}
		
		//! @brief 順伝播計算(ベクトル入出力訓練版)
		//! @param[in]	zprev	前の層からの入力ベクトル
		//! @param[in]	z		出力ベクトル
		void CalcForwardForTraining(const Matrix<1,N>& zprev, Matrix<1,P>& z){
			u = W*zprev + b;						// 重み乗算加算とバイアス加算
			ActivationFunctions::f<AF,1,P>(u, z);	// 活性化関数
			
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトする場合
				z = z & DropMask;	// 該当ノードをドロップアウト
			}
		}
		
		//! @brief 順伝播計算(ベクトル入出力推定版)
		//! @param[in]	zprev	前の層からの入力ベクトル
		//! @param[in]	z		出力ベクトル
		void CalcForwardForEstimation(const Matrix<1,N>& zprev, Matrix<1,P>& z){
			u = W*zprev + b;						// 重み乗算加算とバイアス加算
			ActivationFunctions::f<AF,1,P>(u, z);	// 活性化関数
			
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトしていて推定計算する場合
				z = z*DropRate;		// ドロップアウト率で補正
			}
		}
		
		//! @brief 順伝播計算(ミニバッチ訓練版)
		//! @param[in]	Zprev	前の層からの入力行列
		//! @param[in]	Z		出力行列
		void CalcForwardForTraining(const Matrix<M,N>& Zprev, Matrix<M,P>& Z){
			U = W*Zprev + b*lT;						// 重み乗算加算とバイアス加算
			ActivationFunctions::f<AF,M,P>(U, Z);	// 活性化関数
			
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトする場合
				Z = Z & (DropMask*lT);	// 該当ノードをドロップアウト
			}
		}
		
		//! @brief 順伝播計算(ミニバッチ推定版)
		//! @param[in]	Zprev	前の層からの入力行列
		//! @param[in]	Z		出力行列
		void CalcForwardForEstimation(const Matrix<M,N>& Zprev, Matrix<M,P>& Z){
			U = W*Zprev + b*lT;						// 重み乗算加算とバイアス加算
			ActivationFunctions::f<AF,M,P>(U, Z);	// 活性化関数
			
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトしていて推定計算する場合
				Z = Z*DropRate;		// ドロップアウト率で補正
			}
		}
		
		//! @brief 入力層と内部層(隠れ層)用の誤差行列の計算
		//! @param[in]	WDeltaNext	後ろ側の層からの重み誤差行列
		//! @param[out]	WDelta		重み誤差行列
		void CalcDelta(const Matrix<M,P>& WDeltaNext, Matrix<M,N>& WDelta){
			ActivationFunctions::fp<AF,M,P>(U, fpU);	// 活性化関数の微分を通す計算
			Delta = fpU & WDeltaNext;				// アダマール積の計算
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトする場合
				Delta = Delta & (DropMask*lT);		// 該当ノードをドロップアウト
			}
			WDelta = tp(W)*Delta;					// 前の層に渡すための重み誤差計算
		}
		
		//! @brief 出力層用の誤差行列の計算
		//! @param[in]	Y		出力値行列
		//! @param[in]	D		目標値行列
		//! @param[out]	WDelta	重み誤差行列
		void CalcDeltaForOutputLayer(const Matrix<M,P>& Y, const Matrix<M,P>& D, Matrix<M,N>& WDelta){
			Delta = Y - D;			// 誤差計算
			WDelta = tp(W)*Delta;	// 前の層に渡すための重み誤差計算
		}
		
		//! @brief ドロップアウトマスクの計算
		void CalcDropout(void){
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトするときのみ下記を計算
				DropRand.GetRandomMatrix(DropMask);	// 乱数生成
				
				// ドロップアウト率よりも大きければマスクを0にする
				for(size_t i = 0; i < P; ++i){
					if(DropRate < DropMask.GetElement(1,i+1)){
						DropMask.SetElement(1,i+1, 0);
					}else{
						DropMask.SetElement(1,i+1, 1);
					}
				}
			}
		}

		//! @brief 重み行列とバイアスベクトルの更新
		//! @param[in]	Zprev	前の層からの入力行列
		void UpdateWeight(const Matrix<M,N>& Zprev){
			dW = 1.0/M*Delta*tp(Zprev);	// 重み更新差分値
			db = 1.0/M*Delta*l;			// バイアス更新差分値
			
			// 勾配降下法の指定に従ってコンパイル時条件分岐
			if constexpr(GD == NnDescentTypes::SGD)		CalcSGD(dW, db);		// 確率的勾配降下法
			if constexpr(GD == NnDescentTypes::MOMENTUM)CalcMomentumSGD(dW, db);// モーメンタム確率的勾配降下法
			if constexpr(GD == NnDescentTypes::ADAGRAD)	CalcAdaGrad(dW, db);	// AdaGrad勾配降下法
			if constexpr(GD == NnDescentTypes::RMSPROP)	CalcRMSprop(dW, db);	// RMSprop勾配降下法
			if constexpr(GD == NnDescentTypes::ADADELTA)CalcAdaDelta(dW, db);	// AdaDelta勾配降下法
			if constexpr(GD == NnDescentTypes::ADAM)	CalcAdam(dW, db);		// Adam勾配降下法
		}
		
		//! @brief 重み行列の表示
		void DispWeight(void){
			PrintMatrix(W, "% 16.14e");
		}
		
		//! @brief バイアスベクトルの表示
		void DispBias(void){
			PrintMatrix(b, "% 16.14e");
		}
		
		//! @brief パーセプトロン設定値の表示
		void DispSettings(void){
			printf("xbar     = % 16.14e\n", xbar);		// データセット正規化(標準化)用の平均値
			printf("sigma    = % 16.14e\n", sigma);		// データセット正規化(標準化)用の標準偏差
			printf("eps      = % 16.14e\n", eps);		// 更新ゲイン(SGD, Momentum, AdaGrad, RMSprop)
			printf("alph     = % 16.14e\n", alph);		// 更新ゲイン(Momentum, RMSprop, AdaDelta)
			printf("bet      = % 16.14e\n", bet);		// 更新ゲイン(Adam)
			printf("NearZero = % 16.14e\n", NearZero);	// ゼロ割回避用のゼロに近い値
			printf("DropRate = % 16.14e\n", DropRate);	// ドロップアウト率
		}
		
		//! @brief 訓練/テスト誤差を返す関数
		//! @param[in]	Y	出力行列
		//! @param[in]	D	目標値行列
		//! @return	訓練/テスト誤差(Loss)
		double GetLoss(const Matrix<M,P>& Y, const Matrix<M,P>& D){
			// 活性化関数によって損失関数を変える
			double E = 0;
			if constexpr(AF == ActvFunc::SOFTMAX){
				E = GetCrossEntropy(Y, D);						// 損失関数 クロスエントロピー版
			}else{
				Matrix<1,P> DeltaBar;
				Statistics::MeanRow(Y - D, DeltaBar);			// ミニバッチ内の平均ベクトルを計算
				const Matrix<1,1> e = tp(DeltaBar)*DeltaBar/2.0;// 損失関数 内積版
				E = e[1];
			}
			return E;
		}
		
		//! @brief クロスエントロピーを返す関数
		//! @param[in]	Y	出力行列
		//! @param[in]	D	目標値行列
		//! @return	クロスエントロピー誤差
		double GetCrossEntropy(const Matrix<M,P>& Y, const Matrix<M,P>& D){
			return -sumall(D & loge(Y))/(double)M;		// Softmaxのクロスエントロピー
		}
		
		//! @brief 生計測データセットの正規化(標準化)
		//! @param[in]	x	生計測データ
		template <size_t NN, size_t MM>
		void NomalizeDataset(Matrix<NN,MM>& x){
			xbar  = Statistics::Mean(x);				// 平均の計算
			sigma = Statistics::StandardDeviation(x);	// 標準偏差の計算
			x = (x - xbar)/sigma;	// データセットを正規化(標準化)
		}
		
		//! @brief 入力データの正規化(標準化)
		//! @param[in]	x	入力データ
		void NormalizeInput(Matrix<1,P>& x){
			x = (x - xbar)/sigma;	// データセットのときの値に基づいて入力を正規化(標準化)
		}
		
		//! @brief 重み行列とバイアスベクトルをCSVファイルとして出力する関数
		//! @param[in]	WeightName	重み行列CSVファイル名
		//! @param[in]	BiasName	バイアスベクトルCSVファイル名
		void SaveWeightAndBias(const std::string& WeightName, const std::string& BiasName){
			CsvManipulator::SaveFile(W, WeightName);	// 重み行列を名前を付けて保存
			CsvManipulator::SaveFile(b, BiasName);		// バイアスベクトルを名前を付けて保存
		}
		
		//! @brief 重み行列とバイアスベクトルをCSVファイルとして入力する関数
		//! @param[in]	WeightName	重み行列CSVファイル名
		//! @param[in]	BiasName	バイアスベクトルCSVファイル名
		void LoadWeightAndBias(const std::string& WeightName, const std::string& BiasName){
			CsvManipulator::LoadFile(W, WeightName);	// 重み行列を読み込み
			CsvManipulator::LoadFile(b, BiasName);		// バイアスベクトルを読み込み
		}
		
		//! @brief パーセプトロンの設定をCSVファイルに保存する関数
		//! @param[in]	SettingName	設定値CSVファイル名
		void SaveSettings(const std::string& SettingName){
			Matrix<1,7> S;	// 設定格納用ベクトル
			S.Set(xbar, sigma, eps, alph, bet, NearZero, DropRate);	// 設定値を格納
			CsvManipulator::SaveFile(S, SettingName);				// CSVファイルに保存
		}
		
		//! @brief CSVファイルからパーセプトロンの設定を読み込む関数
		//! @param[in]	SettingName	設定値CSVファイル名
		void LoadSettings(const std::string& SettingName){
			Matrix<1,7> S;	// 設定格納用ベクトル
			CsvManipulator::LoadFile(S, SettingName);				// CSVファイルから読み込み
			S.Get(xbar, sigma, eps, alph, bet, NearZero, DropRate);	// 設定値を読み込み
		}
		
	private:
		SingleLayerPerceptron(const SingleLayerPerceptron&) = delete;					//!< コピーコンストラクタ使用禁止
		const SingleLayerPerceptron& operator=(const SingleLayerPerceptron&) = delete;	//!< 代入演算子使用禁止
	
	protected:
		Matrix<1,P> u;		//!< 状態ベクトル
		Matrix<M,P> U;		//!< 状態行列
		Matrix<N,P> W;		//!< 重み行列
		Matrix<1,P> b;		//!< バイアスベクトル
		Matrix<1,M> l;		//!< 1ベクトル
		Matrix<M,1> lT;		//!< 1ベクトルの転置
		Matrix<1,P> delta;	//!< 誤差ベクトル
		Matrix<M,P> Delta;	//!< 誤差行列
		Matrix<1,P> fpu;	//!< 活性化関数の微分を通した後のベクトル
		Matrix<M,P> fpU;	//!< 活性化関数の微分を通した後の行列
		Matrix<N,P> dW;		//!< 重み更新差分値
		Matrix<1,P> db;		//!< バイアス更新差分値
		Matrix<N,P> DW;		//!< 更新ゲイン乗算後の重み更新差分値
		Matrix<1,P> Db;		//!< 更新ゲイン乗算後のバイアス更新差分値
		Matrix<N,P> HW;		//!< AdaGrad, RMSprop, Adam用
		Matrix<1,P>	Hb;		//!< AdaGrad, RMSprop, Adam用
		Matrix<N,P> GW;		//!< AdaDelta, Adam用
		Matrix<1,P>	Gb;		//!< AdaDelta, Adam用
		Matrix<N,P> HWhat;	//!< Adam用
		Matrix<N,P>	GWhat;	//!< Adam用
		Matrix<1,P> Hbhat;	//!< Adam用
		Matrix<1,P>	Gbhat;	//!< Adam用
		size_t SGDcount;	//!< Adam勾配降下法用のカウンタ
		RandomGenerator DropRand;	//!< ドロップアウト用メルセンヌ・ツイスタ
		Matrix<1,P> DropMask;		//!< ドロップアウト用マスクベクトル
		double xbar;		//!< データセット正規化(標準化)用の平均値
		double sigma;		//!< データセット正規化(標準化)用の標準偏差
		double eps;			//!< 更新ゲイン(SGD, Momentum, AdaGrad, RMSprop)
		double alph;		//!< 更新ゲイン(Momentum, RMSprop, AdaDelta)
		double bet;			//!< 更新ゲイン(Adam)
		double NearZero;	//!< ゼロ割回避用のゼロに近い値
		double DropRate;	//!< ドロップアウト率
		
		//! @brief ヴァニラ確率的勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcSGD(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			DW = -eps*DiffW;	// 更新ゲイン乗算後の重み更新差分値
			Db = -eps*Diffb;	// 更新ゲイン乗算後のバイアス更新差分値
			W += DW;			// 重み行列の更新
			b += Db;			// バイアスベクトルの更新
		}
		
		//! @brief モーメンタム確率的勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcMomentumSGD(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			DW = alph*DW - eps*DiffW;	// 更新ゲイン乗算後の重み更新差分値
			Db = alph*Db - eps*Diffb;	// 更新ゲイン乗算後のバイアス更新差分値
			W += DW;					// 重み行列の更新
			b += Db;					// バイアスベクトルの更新
		}
		
		//! @brief AdaGrad勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcAdaGrad(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			HW += DiffW & DiffW;
			Hb += Diffb & Diffb;
			DW = -eps*(DiffW % (sqrte(HW) + NearZero));
			Db = -eps*(Diffb % (sqrte(Hb) + NearZero));
			W += DW;	// 重み行列の更新
			b += Db;	// バイアスベクトルの更新
		}
		
		//! @brief RMSprop勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcRMSprop(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			HW = alph*HW + (1.0 - alph)*(DiffW & DiffW);
			Hb = alph*Hb + (1.0 - alph)*(Diffb & Diffb);
			DW = -eps*(DiffW % (sqrte(HW) + NearZero));
			Db = -eps*(Diffb % (sqrte(Hb) + NearZero));
			W += DW;	// 重み行列の更新
			b += Db;	// バイアスベクトルの更新
		}
		
		//! @brief AdaDelta勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcAdaDelta(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			// 重みのAdaDeltaの計算
			HW = alph*HW + (1.0 - alph)*(DiffW & DiffW);
			DW = -(sqrte(GW + NearZero) % sqrte(HW + NearZero)) & DiffW;
			GW = alph*GW + (1.0 - alph)*(DW & DW);
			
			// バイアスのAdaDeltaの計算
			Hb = alph*Hb + (1.0 - alph)*(Diffb & Diffb);
			Db = -(sqrte(Gb + NearZero) % sqrte(Hb + NearZero)) & Diffb;
			Gb = alph*Gb + (1.0 - alph)*(Db & Db);
			
			W += DW;	// 重み行列の更新
			b += Db;	// バイアスベクトルの更新
		}
		
		//! @brief Adam勾配降下法
		//! @param[in]	DiffW	重み更新差分値
		//! @param[in]	Diffb	バイアス更新差分値
		void CalcAdam(const Matrix<N,P>& DiffW, const Matrix<1,P>& Diffb){
			++SGDcount;	// イタレーションカウンタ
			
			// 重みのAdamの計算
			HW = alph*HW + (1.0 - alph)*DiffW;
			GW = bet*GW + (1.0 - bet)*(DiffW & DiffW);
			HWhat = HW/(1.0 - pow(alph, SGDcount));
			GWhat = GW/(1.0 - pow(bet,  SGDcount));
			DW = -eps*( HW % (sqrte(GWhat) + NearZero) );
			
			// バイアスのAdamの計算
			Hb = alph*Hb + (1.0 - alph)*Diffb;
			Gb = bet*Gb + (1.0 - bet)*(Diffb & Diffb);
			Hbhat = Hb/(1.0 - pow(alph, SGDcount));
			Gbhat = Gb/(1.0 - pow(bet,  SGDcount));
			Db = -eps*( Hb % (sqrte(Gbhat) + NearZero) );
			
			W += DW;	// 重み行列の更新
			b += Db;	// バイアスベクトルの更新
		}
};
}

#endif

