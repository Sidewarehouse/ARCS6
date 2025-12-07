% ARCS-Matrix行列演算チェック用MATLABコード
% Yokokura, Yuki 2025/01/17
clc;
clear;
format short
%format longE

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
disp '★ ソート関連の関数'
A = [10 -12 4 8; 6 -9 8 0; 2 3 11 -2; 1 1 9 3];
sortA1 = sort(A)

fprintf('\n');
disp '★ 要素ごとの数学関数'
Ax1 = [
		     1.0,      2.0,      3.0,    4.0 ;
		    pi/4,     pi/2, 3.0*pi/4,     pi ;
		5.0*pi/4, 3.0*pi/2, 7.0*pi/4, 2.0*pi ]
Yx1 = exp(Ax1)
Yx6 = tan(Ax1)
Yx7 = tanh(Ax1)

fprintf('\n');
disp '★ 行列間のコピー操作関連の関数'
Fx = [
		10, 20, 30, 40, 50, 60;
		70, 80, 90, 10, 11, 12;
		13, 14, 15, 16, 17, 18;
		19, 20, 21, 22, 23, 24;
		25, 26, 27, 28, 29, 30;
		31, 32, 33, 34, 35, 36;
		37, 38, 39, 40, 41, 42
]
Acpy1(1:10,1:10) = 0;
Acpy1(5:9, 7:8) = Fx(2:6, 3:4)

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
disp '★ 行列式det関連の関数'
Adet1 = [1 -2 4; -5 2 0; 1 0 3]
det(Adet1)
det(A)
det(Ax)
det(Acomp1)

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
[Qqr3_, Rqr3_] = qr(Acomp1)
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
Qqr3 = Qqr3*Di;
Rqr3 = D*Rqr3;

%ArcsMatrixでの複素数QR分解の結果(横長行列)
Acmpx2
[Qqr4_, Rqr4_] = qr(Acmpx2)
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
Qqr4 = Qqr4*Di;
Rqr4 = D*Rqr4;

%ArcsMatrixでの複素数QR分解の結果(縦長行列)
Acmpx2t = Acmpx2'
[Qqr5_, Rqr5_] = qr(Acmpx2t)
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

% ↓実装確認用の残骸↓
Aqr9 = [
	0, 0, 0 ;
	1, 0, 0 ;
	0, 1, 0 ];
A = Acomp1
[Q,R] = qr(A)
[m,n]=size(A);
shouldpermut=false;
isrealA = isreal(A);
notreal = ~all(isrealA(:));
Q = eye(m);
R=A;
P=eye(n);
cols=1:n;
numloop=(m*(m<n)+n*(m>=n));
for k=1:numloop
	col=k;
	if shouldpermut
		%nor2 = sumsq(R(k:m,k:n),1);
		nor2 = sum( R(k:m,k:n).^2 );
		[~,kk] = max(nor2);
		kn = (k:n);
		col = kn(kk);
	end
	a = R(k:m,col);
	if notreal
		alpha = -exp( angle(a(1))*1i )*sqrt(a'*a);
	else
		alpha = sqrt(a'*a);
	end
	e = eye(m-k+1,1)*alpha;
	u = a - e;
	sumsqu = u'*u;
	if sumsqu~=0
		v = u/sumsqu^.5;
	else
		v = u;
	end
	Qk=[eye(k-1) zeros(k-1,m-k+1);zeros(m-k+1,k-1) eye(m-k+1) - 2*v*v'];
	R=Qk*R;
	if col~=k
		Rt=R(:,k);
		R(:,k)=R(:,col);
		R(:,col)=Rt;
		ct=cols(k);
		cols(k)=cols(col);
		cols(col)=ct;
	end
	Q = Q*Qk';
end
Q
R
P = P(:,cols);

%return;	% スクリプトをココで中断

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
[Us5_, Ss5_, Vs5_] = svd(Acmpx2)
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
% しかし、不明な点が多々あり、今後の課題
% A = U S V' = U D2 Di2 S (V D3 Di3)' = U D2 Di2 S Di3' (V D3)' = U D2 S (V D3)' = Un S Vn'
% Di2*S*Di3'    = S
% D2*Di2*S*Di3' = D2*S
%        S*Di3' = D2*S
%          Di3' = S\D2*S
%          Di3  = (S\D2*S)'
%{
d = 1./sign(Us5(1,:))
D2  = diag(d)
Di2 = diag(1./d)
x = (0.9075 + 0.4200i)	% 謎の係数
x = Vs5_(:,3)./Vs5(:,3)	% 謎の係数
Di3 = (Ss5\D2*Ss5)'
Di3(3,3) = 1
D3 = diag(1./diag(Di3))
Z = Ss5*Di3' - D2*Ss5
I = Di3*D3
V  = Vs5
VD = Vs5*D3
det(V)
1/det(VD)	% 謎の係数が計算できたとおもいきや微妙に異なる、不明
%}

%ArcsMatrixでの特異値分解の結果(長い横長行列)
Acomp2 = [
	7.0000 + 8.0000i,   5.0000 + 2.0000i,   7.0000 + 5.0000i,   1.0000 + 6.0000i,   1.0000 + 5.0000i ;
	5.0000 + 7.0000i,   1.0000 + 6.0000i,   1.0000 + 9.0000i,   5.0000 + 8.0000i,   8.0000 + 4.0000i ]
[Us6_, Ss6_, Vs6_] = svd(Acomp2)
Us6 = [
  -0.4248 -  0.4855i,  -0.3786 +  0.6638i  ;
  -0.4101 -  0.6448i,   0.4009 -  0.5053i  ]
Ss6 = [
  23.8292 -  0.0000i,   0.0000 +  0.0000i,  -0.0000 -  0.0000i,  -0.0000 +  0.0000i,  -0.0000 -  0.0000i  ;
   0.0000 +  0.0000i,   8.5539 +  0.0000i,  -0.0000 -  0.0000i,   0.0000 +  0.0000i,   0.0000 +  0.0000i  ]
Vs6 = [
  -0.5632 -  0.0148i,   0.1318 +  0.2738i,  -0.3785 +  0.1004i,  -0.3595 -  0.2983i,  -0.4022 -  0.2386i  ;
  -0.3094 +  0.0100i,  -0.3737 +  0.1362i,  -0.5533 -  0.0898i,   0.0272 +  0.3703i,   0.3643 +  0.4015i  ;
  -0.4874 +  0.0744i,  -0.4066 +  0.2836i,   0.7090 +  0.0034i,  -0.0427 +  0.0353i,   0.0612 +  0.0405i  ;
  -0.4426 +  0.0890i,   0.1830 -  0.3272i,  -0.0141 -  0.0894i,   0.7456 -  0.0144i,  -0.2795 +  0.1150i  ;
  -0.3656 -  0.0789i,   0.4824 -  0.3612i,   0.1153 -  0.0912i,  -0.2613 -  0.1280i,   0.6251 -  0.0012i  ]
% MATLABに準拠させるための等価変換
% しかし、不明な点が多々あり、今後の課題
%{
x = Vs6_(:,3)./Vs6(:,3)	% 謎の係数
x = Vs6_(:,4)./Vs6(:,4)	% 謎の係数
x = Vs6_(:,5)./Vs6(:,5)	% 謎の係数
%}

fprintf('\n');
disp '★ 階数関連の関数'
Arank1 = [3 2 4; -1 1 2; 9 5 10]
rankA1 = rank(Arank1)
Arank2 = [10 0 0 0; 0 25 0 0; 0 0 34 0; 0 0 0 1e-15]
rankA2 = rank(Arank2)

fprintf('\n');
disp '★ コレスキー分解関連の関数'
Achol1 = [
	 7, -2,  9 ;
	-2,  6, -3 ;
	 9, -3,  3
]
[Lldl1, Dldl1] = ldl(Achol1)
Acomp3 = [
	4.0 + 6.0i,  1.0 - 3.0i,  5.0 + 2.0i ;
	1.0 - 3.0i, -7.0 - 6.0i,  7.0 - 1.0i ;
	5.0 + 2.0i,  7.0 - 1.0i, -5.0 - 3.0i ]
%[Lldl2, Dldl2] = ldl(Acomp3)	% 複素数のLDL分解はMATLABでは非対応
%chol(Achol1)					% 対称正定値行列ではないのでエラー
Achol2 = [
    2.3370   -0.0127   -0.6299 ;
   -0.0127    2.3370   -0.0127 ;
   -0.6299   -0.0127    2.3370 ]
Rchol1 = chol(Achol2)

fprintf('\n');
disp '★ 線形方程式関連の関数'
Aslv1 = [
	1,  1,  1 ;
	2,  3, -2 ;
	3, -1,  1 ]
bslv1 = [
	9 ;
	5 ;
	7 ]
xslv1 = linsolve(Aslv1, bslv1)
Bslv2 = [
	9,  2,  7 ;
	5,  6, -7 ;
	7, -3,  5 ]
Xslv2 = linsolve(Aslv1, Bslv2)
Aslv3 = [
	1,  1,  1 ;
	2,  3, -2 ;
	3, -1,  1 ;
	7,  5,  4 ]
bslv3 = [
	9 ;
	5 ;
	7 ;
	3 ]
xslv3 = linsolve(Aslv3, bslv3)
Bslv4 = [
	9,  2,  7 ;
	5,  6, -7 ;
	7, -3,  5 ;
	3,  1,  2 ]
Xslv4 = linsolve(Aslv3, Bslv4)
xslv5 = linsolve(Aslv3', bslv1)		% 縦長の場合、ArcsMatrixとは異なる解が得られる
Xslv6 = linsolve(Aslv3', Bslv4')	% 縦長の場合、ArcsMatrixとは異なる解が得られる

fprintf('\n');
disp '★ 逆行列関連の関数'
Ainv1 = [1 0 2; -1 5 0; 0 3 -9]
Yinv1 = inv(Ainv1)
Yinv3 = inv(Acomp1)
Yinv4 = pinv(Aslv3)
Yinv5 = pinv(Aslv3')

fprintf('\n');
disp '★ Hessenberg分解関連の関数'
Asch1 = [
 -149    -50   -154 ;
  537    180    546 ;
  -27     -9    -25 ]
[Phes1, Hhes1] = hess(Asch1)
[Phes2, Hhes2] = hess(Aqr1)
[Phes3, Hhes3] = hess(Acomp1)
% ↓実装確認用の残骸↓
%{
%A = Asch1;
%A = Aqr1
A = Acomp1
%A = rand(4)
%A(2,1) = 0
Aiscomplex = ~isreal(A);
[n,m] = size(A);
I = eye(n);
P = I;
H = A;
u = zeros(n,1)
for k=1:n-2
	k
	r = [H(k+1:n,k); zeros(k,1)]
	if Aiscomplex
		u(1) = -exp(angle(r(1))*1i)*sqrt(r'*r)
	else
		if r(1) == 0
			u(1) = -sqrt(r'*r);
		else
			u(1) = -sign(r(1))*sqrt(r'*r);
		end
	end
	v = r - u
	v = v/sqrt(v'*v);
	vv = v*v';
	I_2vv = I - 2*v*v';
	W = [
		eye(k),       zeros(k,n-k)      ;
		zeros(n-k,k), I_2vv(1:n-k,1:n-k) ]
	H = W*H*W';
	P = P*W';
end
P
H
norm(Acomp1 - P*H*P')
%}

fprintf('\n');
disp '★ Schur分解関連の関数'
[Usch1, Tsch1] = schur(Asch1)
%[Usch2, Tsch2] = schur(Acomp1)
%[Usch3, Tsch3] = schur(magic(4))
%UTU = Usch3*Tsch3*Usch3'
Asch2 = [
	 4, -8,  1;
	 1, -5,  1;
	-9,  8, -6
]
[Usch2, Tsch2] = schur(Asch2)
Asch3 = [
	 1,  2,  3 ;
	-5,  9, -1 ;
	 2,  6,  8
]
[Usch2, Tsch2] = schur(Asch3)
% ↓実装確認用の残骸↓
%{
A = Asch3
[n, m] = size(A);
[P, H] = hess(A);
k = n;
U = P
S = H
tol = 1e-14;
while k > 1
	k = nnz(diag(S,-1)) + 1
	if k > 1
		a = ( S(k-1,k-1)+S(k,k) + sqrt( (S(k-1,k-1) + S(k,k))^2 - 4*(S(k-1,k-1)*S(k,k) - S(k-1,k)*S(k,k-1)) ) )/2
		disp 'a*I = '
		a*eye(k)
		W = S(1:k,1:k) - a*eye(k)
		[Q R] = qr(W);
		Q
		V = [               Q, zeros(k-1+1,n-k) ;
		     zeros(n-k,k-1+1),         eye(n-k) ]
		U = U*V
		S = V'*S*V;
		Z = 1 - tril(abs(S) < tol, -1)
		S = S.*Z
	end
end
U
S
U*S*U'
%}

fprintf('\n');
disp '★ 固有値関連の関数'
Aeig1 = gallery("lehmer",4)
veig1 = eig(Aeig1)
Aeig2 = gallery("circul",3)
veig2 = eig(Aeig2)
Aeig3 = [3 1 0; 0 3 1; 0 0 3]
veig3 = eig(Aeig3)
Aeig4 = [-3, -4,  2; -7,  1, -5; 6, -7,  3]
veig4 = eig(Aeig4)
Aeig5 = [10,  -8,   5; -8,   9,   6; -1, -10,   7]
veig5 = eig(Aeig5)
[Veig1, Deig1] = eig(Aeig1)
[Veig2, Deig2] = eig(Aeig2)
[Veig2, Deig2] = eig(Aeig3)

% ↓実装確認用の残骸↓
%{
A = Aeig3;
[m, n]=size(A);
I = eye(n);
[U, S] = schur(A,'complex');
l = diag(S);
V = zeros(n);
for k = 1:n
	(A - l(k)*I)'
	[Q,~] = qr( (A - l(k)*I)' );
	Q
	V(:,k) = Q(:,n);
end
V
D = diag(l)
ll_AV_VD_ll = norm(A*V - V*D)
%}

fprintf('\n');
disp '★ その他の「積」関連の関数'
Lkrn1 = [1,2;3,4];
Rkrn1 = [0,5;6,7];
Ykrn1 =  kron(Lkrn1, Rkrn1)	% クロネッカー積
Lkrn2 = [
	1,  1,  1 ;
	2,  3, -2 ;
	3, -1,  1
];
Rkrn2 = [
	1, 2;
	3, 4;
	5, 6;
	7, 8
];
Ykrn2 = kron(Lkrn2, Rkrn2)
Ykrn2 = kron(Rkrn2, Lkrn2)

fprintf('\n');
disp '★ 行列指数関数'
Aexp1 = [
	1,  1,  0 ;
	0,  0,  2 ;
	0,  0, -1 ];
Yexp1 = expm(Aexp1)
format longE
Yexp1 = expm(Aexp1)
format short
Yexp2x = expm(Lkrn2)

fprintf('\n');
disp '★ 統計処理関数'
Ast1 = [
		0, 1, 1;
		2, 3, 2;
		1, 3, 2;
		4, 2, 2];
yst1 = mean(Ast1)
yst2 = mean(Ast1, 2)
meanAst1 = mean(mean(Ast1))
Ast2 = [4 -2 1; 9 5 7];
yst3 = var(Ast2)
yst4 = var(Ast2, 0 , 2)
varyst3 = var(yst3)
varyst4 = var(yst4)
varAst2 = var(Ast2, 0 , 'All')
Ast3 = [6 4 23 -3; 9 -10 4 11; 2 8 -5 1];
yst5 = std(Ast3)
yst6 = std(Ast3, 0, 2)
varyst5 = std(yst5, 0, 'All')
varyst6 = std(yst6, 0, 'All')
varAst3 = std(Ast3, 0, 'All')

fprintf('\n');
disp '★ その他の行列関連関数'
vpol1 = rot90([2.0 + 1.0i; 3.0 + 2.0i; 5.0 + 3.0i], 2)
vcoef1 = poly(vpol1)

