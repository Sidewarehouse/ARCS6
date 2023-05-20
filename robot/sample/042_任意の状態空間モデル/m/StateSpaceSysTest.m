% StateSpaceSystem�N���X�̃e�X�g�p�X�N���v�g
% 2022/02/20 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% �e�X�g�Ώۂ̏�ԋ�ԃ��f��
A = [
	 0,  1 ;
	-5, -2
];
b = [
	0 ;
	3
];
C = [
	1, 0
];
D = [
	1
];
System1 = ss(A,b,C,D)

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
u = CsvData(:,2);
y = CsvData(:,3);
clear CsvData;
tlen = length(t);

% �^�l�̌v�Z
Ts = t(tlen)/(tlen-1);
ti = 0:Ts:t(tlen);
yi = lsim(System1, u, ti);

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
h=stairs(t, u, 'k:');
	set(h,'linewidth',2);
hold on;
h=stairs(t, y, 'r');
	set(h,'linewidth',3);
h=stairs(t, yi, 'k');
	set(h,'linewidth',1);
hold off;
xlabel('Time [s]','FontSize',12);
ylabel('Input and Output [-]','FontSize',12);
set(gca,'FontSize',12);
grid on;
%axis([0 10 -inf inf]);
legend('Input','Output','Ideal','Location','SouthEast','Orientation','Vertical');

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));
