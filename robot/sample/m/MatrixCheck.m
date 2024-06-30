% 行列演算チェック用MATLABコード
% Yokokura, Yuki 2024/06/27
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
fprintf('\n');

disp '★ ノルム関連の関数(行列版)'
Ax1 = [
	     1.0,      2.0,      3.0,    4.0 ;
	    pi/4,     pi/2, 3.0*pi/4,     pi ; 
	5.0*pi/4, 3.0*pi/2, 7.0*pi/4, 2.0*pi 
]
fprintf('norm(Ax1) = %f\n', norm(Ax1));
fprintf('norm(Ax1, 1) = %f\n',norm(Ax1, 1)); 
fprintf('norm(Ax1, Inf) = %f\n', norm(Ax1, Inf));
Acmpx2 = [
  1 + 1j,  3 + 4j,  3 - 4j ;
 -1 + 1j, -3 + 4j, -3 - 4j
]
fprintf('norm(Acmpx2) = %f\n', norm(Acmpx2));
fprintf('norm(Acmpx2, 1) = %f\n', norm(Acmpx2, 1));
fprintf('norm(Acmpx2, Inf) = %f\n\n', norm(Acmpx2, Inf));

disp '★ ノルム関連の関数(ベクトル版)'
k1 = [
  2.000 ;
  5.000 ;
 -9.000 ;
  4.000 ;
  7.000 ;
 -1.000 ;
  6.000
]
fprintf('norm(k1) = %f\n', norm(k1));
fprintf('norm(k1, 1) = %f\n', norm(k1, 1));
fprintf('norm(k1, Inf) = %f\n', norm(k1, Inf));
k1cmpx = sqrt(-k1)
fprintf('norm(k1cmpx) = %f\n', norm(k1cmpx));
fprintf('norm(k1cmpx, 1) = %f\n', norm(k1cmpx, 1));
fprintf('norm(k1cmpx, Inf) = %f\n\n', norm(k1cmpx, Inf));

disp '★ LU分解関連の関数'
A = [
	10, -7,  0 ;
	-3,  2,  6 ;
	 5, -1,  5
]
[L, U, P] = lu(A)
Ax = [
	 5,  6,  7 ;
	 8,  9,  0 ;
	 1,  2,  3
]
[L, U] = lu(Ax)
Acomp1 = [
	4 + 6i,  1 - 3i,  5 + 2i ;
	8 - 5i, -7 - 6i,  7 - 1i ;
	9 + 9i, -7 - 5i, -5 - 3i
]
[L, U, P] = lu(Acomp1)

fprintf('\n');
disp '★ QR分解関連の関数'
Aqr1 = magic(5)
[Qqr1, Rqr1] = qr(Aqr1)
fprintf('norm(Aqr1 - Qqr1*Rqr1) = %e\n\n', norm(Aqr1 - Qqr1*Rqr1) );
Aqr2 = [
	12, -51,   4, 39;
	 6, 167, -68, 22;
	-4,  24, -41  11;
]
[Qqr2, Rqr2] = qr(Aqr2)
[Qx2, Rx2] = qr(Aqr2')

%ArcsMatrixでの複素数QR分解の結果(正方行列)
Acomp1
[Qqr3, Rqr3] = qr(Acomp1)
Qqr3 = [
  -0.414 -  0.000i,   0.165 -  0.101i,   0.888 +  0.051i ;
  -0.016 +  0.542i,  -0.751 -  0.064i,   0.104 +  0.356i ;
  -0.717 +  0.143i,  -0.001 +  0.628i,  -0.259 -  0.064i
];
Rqr3 = [
  -9.656 - 14.483i,   0.749 +  9.719i,   0.430 -  1.737i ;
   0.000 +  0.000i,   2.984 +  8.067i,  -6.449 +  5.182i ;
   0.000 +  0.000i,   0.000 +  0.000i,   6.401 -  0.623i
];
D  = diag(1./sign(diag(Rqr3)));
Di = diag(sign(diag(Rqr3)));
Q = Qqr3*Di
R = D*Rqr3

%ArcsMatrixでの複素数QR分解の結果(横長行列)
Acmpx2
[Qqr4, Rqr4] = qr(Acmpx2)
Qqr4 = [
  -0.707 +  0.000i,   0.000 -  0.707i ;
   0.000 -  0.707i,  -0.707 -  0.000i 
];
Rqr4 = [
  -1.414 -  1.414i,  -4.950 -  4.950i,   0.707 +  0.707i, ;
   0.000 +  0.000i,  -0.707 -  0.707i,   4.950 +  4.950i
];
D  = diag(1./sign(diag(Rqr4)));
Di = diag(sign(diag(Rqr4)));
Q = Qqr4*Di
R = D*Rqr4

%ArcsMatrixでの複素数QR分解の結果(縦長行列)
Acmpx2.'
[Qqr5, Rqr5] = qr(Acmpx2.')
Qqr5 = [
  -0.196 +  0.000i,  -0.199 +  0.041i,   0.067 +  0.957i ;
  -0.686 +  0.098i,  -0.663 +  0.040i,  -0.020 -  0.279i ;
   0.098 -  0.686i,  -0.144 +  0.705i,  -0.003 -  0.040i 
];
Rqr5 = [
  -5.099 +  5.099i,  -1.177 +  1.569i ;
  -0.000 -  0.000i,   5.239 +  4.551i ;
  -0.000 -  0.000i,   0.000 -  0.000i 
];
D  = diag(1./sign(diag(Rqr5)))
Di = diag(sign(diag(Rqr5)))
%{ % Qが補正できないので今後の課題
D = [D ; [0, 0]];
D = [D [0 ; 0; 1]]
Di = [Di ; [0, 0]];
Di = [Di [0 ; 0; -1]]
Q = Qqr5*Di
R = D*Rqr5
%}
%A = QR = ((QR)t)t = ( ((QDi)(DR))t )t = ( (DR)t (QDi)t )t = ( RtDt DitQt )t

%{
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

