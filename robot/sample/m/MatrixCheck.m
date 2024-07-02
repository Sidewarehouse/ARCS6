% 行列演算チェック用MATLABコード
% Yokokura, Yuki 2024/07/02
clc;
clear;
format short

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
% MATLABに準拠させるための等価変換
D  = diag(1./sign(diag(Rqr3)));
Di = diag(sign(diag(Rqr3)));
Q = Qqr3*Di;
R = D*Rqr3;


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
% MATLABに準拠させるための等価変換
D  = diag(1./sign(diag(Rqr4)));
Di = diag(sign(diag(Rqr4)));
Q = Qqr4*Di;
R = D*Rqr4;

%ArcsMatrixでの複素数QR分解の結果(縦長行列)
Acmpx2t = Acmpx2'
[Qqr5, Rqr5] = qr(Acmpx2t)
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
% MATLABに準拠させるための等価変換
% しかし、Qが補正できないので今後の課題
%{
D  = diag(1./sign(diag(Rqr5)));
Di = diag(sign(diag(Rqr5)));
D = [D ; [0, 0]];
D = [D [0 ; 0; 1]]
Di = [Di ; [0, 0]];
Di = [Di [0 ; 0; -1]]
Q = Qqr5*Di
R = D*Rqr5
%}

fprintf('\n');
disp '★ SVD特異値分解関連の関数'
As1 = [
	1, 2 ;
	3, 4 ;
	5, 6 ;
	7, 8
]
[Us1, Ss1, Vs1] = svd(As1)
[Us2, Ss2, Vs2] = svd(As1')

% SVDのQR分解部分のデバッグ(Aが横長の場合)
%{
U = eye(2);
V = eye(4);
S = As1';
for l = 1:5
	[Qn, Snm] = qr(S');
	S = Snm';
	[Qm, S] = qr(S);
	
	U = U*Qm;
	V = V*Qn;
end
%}

As2 = [
	 2,  0,  2 ;
	 0,  1,  0 ;
	 0,  0,  0
]
[Us3, Ss3, Vs3] = svd(As2)

As2 = [
	 1,  1,  3,
	-5,  6, -3, 
	 7, -2,  9
]
[Us3, Ss3, Vs3] = svd(As2)

%ArcsMatrixでの特異値分解の結果(正方行列)
[Us4, Ss4, Vs4] = svd(Acomp1)
Us4 = [
   0.1996 +  0.3438i,   0.2459 -  0.1672i,   0.8580 -  0.1319i ;
   0.6064 -  0.1480i,   0.1219 -  0.6034i,  -0.1649 +  0.4519i ;
   0.5746 +  0.3495i,  -0.5565 +  0.4721i,  -0.0034 +  0.1226i ];
Ss4 = [
  20.5503 +  0.0000i,   0.0000 -  0.0000i,  -0.0000 -  0.0000i ;
  -0.0000 +  0.0000i,  12.1589 +  0.0000i,   0.0000 +  0.0000i ;
   0.0000 -  0.0000i,   0.0000 -  0.0000i,   3.8534 +  0.0000i ];
Vs4 = [
   0.8160 -  0.0000i,   0.2642 +  0.2381i,   0.0351 -  0.4543i ;
  -0.4846 +  0.2941i,   0.4154 -  0.0462i,  -0.2318 -  0.6710i ;
   0.1049 +  0.0422i,   0.3058 -  0.7780i,   0.5370 -  0.0000i ];
% MATLABに準拠させるための等価変換
% A = U S V' = U D Di S (V D Di)' = U D Di S Di' (V D)' = U D S (V D)' = Un S Vn'
% Di = diag(1./d);
% D*Di' = I
d = 1./sign(Vs4(1,1:3));
D  = diag(d);
Us4 = Us4*D
Vs4 = Vs4*D
Z = Acomp1 - Us4*Ss4*Vs4'

%ArcsMatrixでの特異値分解の結果(横長行列)
[Us5, Ss5, Vs5] = svd(Acmpx2)
Us5 = [
  -0.5000 -  0.5000i,  -0.4629 +  0.5346i ;
  -0.4243 -  0.5657i,   0.5338 -  0.4637i ];
Ss5 = [
   8.1328 -  0.0000i,   0.0000 +  0.0000i,   0.0000 +  0.0000i ;
   0.0000 -  0.0000i,   6.1529 -  0.0000i,   0.0000 +  0.0000i ];
Vs5 = [
  -0.1403 +  0.1217i,  -0.1505 +  0.1507i,   0.0672 +  0.9569i ;
  -0.5521 +  0.4788i,  -0.4399 +  0.4406i,  -0.0196 -  0.2791i ;
   0.4962 -  0.4304i,  -0.5320 +  0.5329i,  -0.0028 -  0.0399i ];
% MATLABに準拠させるための等価変換
% A = U S V' = U D2 Di2 S (V D3 Di3)' = U D2 Di2 S Di3' (V D3)' = U D2 S (V D3)' = Un S Vn'
% Di2*S*Di3'    = S
% D2*Di2*S*Di3' = D2*S
%        S*Di3' = D2*S
%          Di3' = S\D2*S
%          Di3  = (S\D2*S)'
norm(Vs5(:,3))
d = 1./sign(Us5(1,:))
D2  = diag(d)
Di2 = diag(1./d)
%Di3 = [ diag(diag(D2*Ss5)./diag(Ss5)), zeros(2,1) ]
Di3 = (Ss5\D2*Ss5)'
Z = Ss5*Di3' - D2*Ss5
%{
D3 = diag(1./diag(Di3))
D3(3,3) = 0
D3*Di3
Vs5*D3
%}
%Us5 = Us5*D

%{
Us5(:,1) = d(1)*Us5(:,1);
Us5(:,2) = d(2)*Us5(:,2);
Us5
Vs5(:,1) = d(1)*Vs5(:,1);
Vs5(:,2) = d(2)*Vs5(:,2);
Vs5(:,3) = (0.9075 + 0.4200i)*Vs5(:,3);	% ←謎のスカラー倍するとMATLABと合う
Vs5
%}

%[Us6, Ss6, Vs6] = svd(Acmpx2')


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

