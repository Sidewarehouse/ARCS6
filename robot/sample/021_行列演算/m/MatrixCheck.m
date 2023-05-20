% �s�񉉎Z�`�F�b�N�pMATLAB�R�[�h
% Yokokura, Yuki 2022/07/27
clc;
clear;

disp '�� �s��錾�ƃZ�b�g�̃e�X�g'
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

disp '�� �c�x�N�g����[]�I�y���[�^�ɂ��v�f�A�N�Z�X�̃e�X�g'
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

disp '�� �s���()�I�y���[�^�ɂ��v�f�A�N�Z�X�̃e�X�g'
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

disp '�� ���Z�q�̃e�X�g'
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


%{
% QR����
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

% Schur����
Asr1 = [
	 4, -8,  1,
	 1, -5,  1,
	-9,  8, -6
]
[Qsr,Usr] = schur(Asr1)
Qsr*Usr*inv(Qsr)

% �N���l�b�J�[�ς̃e�X�g
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

