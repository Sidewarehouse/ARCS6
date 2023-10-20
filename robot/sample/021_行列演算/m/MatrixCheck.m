% 行列演算チェック用MATLABコード
% Yokokura, Yuki 2022/07/27
clc;
clear;

disp '◆ 行列宣言とセットのテスト'
A = [
	1,  1,  1 ;
	2,  3, -2 ;
	3, -1,  1
];
B = [
	5,  1, -3 ;
	3, -1,  7 ;
	4,  2, -5
];
C = A;
size(A)
A
B
C

disp '◆ 縦ベクトルの[]オペレータによる要素アクセスのテスト'
alpha = [
	1 ;
	2 ;
	3 ;
	4 ;
	5
];
size(alpha)
alpha
alpha(3)
alpha(3) = 6;
alpha

disp '◆ 行列の()オペレータによる要素アクセスのテスト'
B(1,1)
B(1,2)
B(1,3)
B(2,1)
B(2,2)
B(2,3)
B(3,1)
B(3,2)
B(3,3)
B(2,3) = 0;
B
B(2,3) = 7;
B

disp '◆ 演算子のテスト'
C = +A
C = -A
C = A + B
C = A + 3
C = A - 3
C = A - B
C = A*B
C = A*2
C = A/2
C = C + A
C = C + 2
C = C - A
C = C - 1
C = C*A
C = C*2
C = A^0
C = A^1
C = A^2
C = A^3
C = A.*B
C = A./B
C = 3 + A
C = 3 - A
C = 2*A

disp '★ ノルム関連の関数(行列版)'
Ax1 = [
	     1.0,      2.0,      3.0,    4.0 ;
	    pi/4,     pi/2, 3.0*pi/4,     pi ; 
	5.0*pi/4, 3.0*pi/2, 7.0*pi/4, 2.0*pi 
]
norm(Ax1)
norm(Ax1, 1)
norm(Ax1, Inf)
Acmpx2 = [
  1 + 1j,  3 + 4j,  3 - 4j ;
 -1 + 1j, -3 + 4j, -3 - 4j
]
norm(Acmpx2)
norm(Acmpx2, 1)
norm(Acmpx2, Inf)

disp '★ ノルム関連の関数(ベクトル版)'
k1 = [
  2.000 ;
  5.000 ;
 -9.000 ;
  4.000 ;
  7.000 ;
 -1.000 ;
  6.000
];
norm(k1)
norm(k1, 1)
norm(k1, Inf)


%{
% QR分解
Aqr = [
	12, -51,   4;
	 6, 167, -68;
	-4,  24, -41;
];
[Qqr, Rqr] = qr(Aqr)

Aqr = [
	12, -51,   4, 39;
	 6, 167, -68, 22;
	-4,  24, -41  11;
];
[Qqr, Rqr] = qr(Aqr)
[Qqr, Rqr] = qr(Aqr.')

Aqr = [
	10, -8,  5;
	-8,  9,  6;
	-1,-10,  7;
];
[Qqr, Rqr] = qr(Aqr)

Acomp3 = [
  4 + 6i,  1 - 3i,  5 + 2*i ;
  8 - 5i, -7 - 6i,  7 - 1*i ;
  9 + 9i, -7 - 5i, -5 - 3*i ;
];
inv(Acomp3)

% Schur分解
Asr1 = [
	 4, -8,  1,
	 1, -5,  1,
	-9,  8, -6
]
[Qsr,Usr] = schur(Asr1)
Qsr*Usr*inv(Qsr)

% クロネッカー積のテスト
kron([1,2;3,4], [0,5;6,7])
Axsvd = [
	1, 2;
	3, 4;
	5, 6;
	7, 8
];
kron(A,Axsvd)
kron(Axsvd,A)
%}

