% �N�[�������C�ƔS�����C�W���̓���X�N���v�g
% 2021/12/10 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% �p�����[�^�ݒ�
Ktn  = 0.31;	% [Nm/A] �g���N�萔
Nave = 6000;	% [-] ���΍����o�̂��߂̈ړ����ϓ_��
wzero = 1e-2;   % [rad/s] ��Ƃ݂Ȃ����x

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
iqref = CsvData(:,2);
wm = CsvData(:,3);
clear CsvData;
tlen = length(t);

% ��푬�x�����o���邽�߂̈ړ����f�B�A���t�B���^
wm_filt = movmedian(wm, [Nave 0]);

% ����Ԃɂ�����e�l�̌��o
j = 0;
for i = 2:tlen
	if iqref(i - 1) ~= iqref(i)			% �d���w�ߒl���ω�������C
		j = j + 1;
		t_stead(j)  = t(i - 1);			% ����Ԃ̎���
		iq_stead(j) = iqref(i - 1);		% ���d��
		wm_stead(j) = wm_filt(i - 1);	% ��푬�x
	end
end
stead_len = length(wm_stead);

% ���xvs�g���N�����ւ̕ϊ�
Ktiq_stead = Ktn*iq_stead;			% [Nm] ����g���N�̌v�Z
Ktiq_stead_abs = abs(Ktiq_stead);	% [Nm] ����]�ƕ���]�̓����𓯈ꎋ����
wm_stead_abs   = abs(wm_stead);		% [rad/s] ����]�ƕ���]�̓����𓯈ꎋ����

% ���x����̏ꍇ�̓t�B�b�e�B���O���f�[�^���珜�O����
% (�ŏ����@�̓��萸�x���グ�邽�߂̑O����)
Ktiq_lsq = Ktiq_stead_abs;
wm_lsq = wm_stead_abs;
for i = stead_len:-1:1,
   if(-wzero < wm_lsq(i) && wm_lsq(i) < wzero)	% ���x����Ȃ�C
       Ktiq_lsq(i) = [];	% �f�[�^�폜
       wm_lsq(i) = [];		% �f�[�^�폜
   end
end
lsq_len = length(wm_lsq);

% �ŏ����@�ɂ�門�C����
A(1:lsq_len,1) = wm_lsq.';
A(1:lsq_len,2) = 1;
b = Ktiq_lsq.';
x = A\b;						% �o�b�N�X���b�V�����Z�q�� A*x = b ��x�ɂ��ĉ���
Dall  = x(1);					% [Nm/(rad/s)] �S�����C�W��(1���֐��̌X��)
Fclmb = x(2);					% [Nm] �N�[�������C��(1���֐��̐ؕ�)
Ktiq_fit = wm_lsq*Dall + Fclmb;	% [Nm] ���C�g���N�t�B�b�e�B���O�����̌v�Z
fprintf('Dall  = %3.2e;	// [Nm/(rad/s)] �S�S�����C�W��\n', Dall);
fprintf('Fclmb = %3.2e;	// [Nm] �N�[�������C��\n', Fclmb);

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','auto');
set(gcf,'color',[1 1 1]);
subplot(2,2,1);
	h=plot(t, iqref, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('q-axis current [A]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(2,2,3);
	h=plot(t, wm, 'k', t, wm_filt, 'r', t_stead, wm_stead, 'go');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Velocity [rad/s]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	legend('Raw','Moving Median','Steady State','Location','NorthWest','Orientation','Vertical');
subplot(2,2,[2,4]);
	h=plot(wm_stead_abs, Ktiq_stead_abs, 'go', wm_lsq, Ktiq_lsq, 'k+', wm_lsq, Ktiq_fit, 'r');
		set(h,'linewidth',2);
	xlabel({'Velocity [rad/s]','(c)'},'FontSize',12);
	ylabel('Motor Torque [Nm]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	legend('Steady State','Non-zero Speed','Fitted by Least Sq.','Location','SouthEast','Orientation','Vertical');
	title( strcat('Dall = ', sprintf('%4.3e',Dall), ' [Nm/(rad/s)]', '   Fclmb = ', sprintf('%4.3e',Fclmb), ' [Nm]') );
saveas(gcf, strcat(FileName,'.png'));

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));
