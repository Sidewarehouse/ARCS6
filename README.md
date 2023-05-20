# ARCS6
ARCS6 - Advanced Robot Control System V6

★ コレは何？ ★
　Linuxで低層レイヤでのロボット制御をするためのフレームワークのようなもの。 
ARCS5.1の後継版。 良くある一般的なフツーのLinuxで，カーネルを一切イジることなく，
さらにユーザ空間上でのコーディングであってもリアルタイム性の確保もしてくれる。
言語はC++17

詳細は → https://www.sidewarehouse.net を参照

★ ARCS6の主な機能 ★
- マルチリアルタイムスレッドの生成・周期実行・破棄
- (RTLinuxやRTAI等々のリアルタイムカーネルを使うことなく，ARCS単体で制御周期100μsのリアルタイム制御を実現！)
- 測定データの記憶と出力(CSV形式及びTAB区切りDAT形式対応)
- アクチュエータパラメータ(位置，電流，推力等)のリアルタイム表示とリアルタイム波形描画
- 任意変数値のリアルタイム表示とリアルタイム波形描画
- オンライン変数書き換え機能
- 一時停止＆再始動機能
- 実際の制御周期の計測と標準偏差，最大値，最小値の計算と表示
- 実際の計算消費時間の計測と表示
- D/Aコンバータ，A/Dコンバータ，エンコーダカウンタ，シリアル通信ボード等々とのインターフェースクラス
- 6軸力覚センサ用シリアル通信インターフェース
- MECHATROLINK III 通信対応(内部版のみ)，EtherCAT通信対応，UNINET通信対応
- デバイスドライバ，カーネルモジュール開発支援機能
- モーションコントロールクラス(積分器，擬似微分器，PI/PD/PID/I-PD等々の制御器，位相補償器，各種オブザーバ，各種フィルタ，レギュレータ等々)
- 様々な機能を持ったクラス(リングバッファ，統計処理，リミッタ，信号発生器，UDP送受信器，メルセンヌ・ツイスタ…等々)
- リアルタイムデバッグプリント機能，イベントログ機能，条件指定緊急停止機能
- 行列＆ベクトルの加減乗演算，冪乗，LU/QR/コレスキー分解，SVD，逆行列，複素数行列，constexpr行列…等々の行列計算
- 連続系状態方程式のA行列とB行列の離散化，任意の連続系伝達関数の応答計算
- 周波数特性測定のためのFRA(Frequency Response Analyzer)クラス
- 単純/単層パーセプトロン，順伝播ニューラルネットワーク

★ ROS/ROS2との違い ★
インバータのゲートを叩くレベルでロボットを動かしている我々パワー研からすると，ROS/ROS2は上位層用である。 
なので，ROSで上位層を，ARCS6で下位層を作るという住み分けが正解。

★ コーディング規約 ★
- 関数名・変数名はアッパーキャメルケース(パスカル)で命名すること。
- 例外的に，数式ライクにした方が見栄えが良い関数名はスネークケースも認める。
- タブインデントは ./lib/ClassBase.hh/cc, ./lib/ClassTemplate.hh/cc, ./lib/FunctionBase.hh/cc を踏襲すること。
- 一時的にコメントアウトしたコードは最終的に消去すること。
- using namespace std; は例外なく使用禁止。
- 定数値は #define でも const ではなく constexpr を使用すること。マクロは使って良い。
- グローバル変数の使用は極力回避すること。
- goto文の使用は極力回避すること。
- new/delete, shared_ptrの使用は極力回避し，unique_ptrをなるべく使用すること。
- Cの古典的ポインタの使用は極力回避し，C++の現代的「参照」をなるべく使用すること。
- コメント文は母国語の「口語的表現」でなるべく付けて，コードとの相違がないようにすること。
- プロトタイプ宣言の引数変数は省略せずに書くこと。
- Efficient C++, More Efficient C++ にできるだけ準拠すること。
- コードを書いたらCppcheck等々により静的解析をすること。
- 物理定数，物理変数はすべてSI単位系で記述すること。

★ ライセンス ★

ARCS6は2023からMIT Licenseにしました。GPLはやめて、BSDからも移行しました。

Copyright (c) 2023 Yokokura, Yuki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


