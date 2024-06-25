% MultiPlot �}���`(�}���^�C)�v���b�g�N���X
% 2024/04/03 Yokokura, Yuki
%
% �{�v���b�g�̓���
% �E�]�������S�ɐ��䉺�ɒu���鑽�i�v���b�g
% �@(�f�t�H��subplot�ł͗]����肷���Ę_���p�̐}�ɂ��)
% �E�t�H���g���E�T�C�Y�A���̐F�E�����A�ȂǁA���܂�؂����ݒ�����Ȃ��ėǂ�
%
% �@ �T�C�Y��`�p�����[�^�̈ʒu�֌W�́��̒ʂ�
%{
    0%                                                  100%
100% +--------------------- FigWidth ----------------------+
     |       |TopMargin
     |       +        +--------------------------+
     |       |        |       PlaneNum = 2       |
     |       |        |PlaneHeight               |
     |       |        |                          |
     |       |        +--------------------------+
     |       |        +--------------------------+ <- XaxisSubMargin
     |Fig    |All     |       PlaneNum = 1       |
     |Height |Planes  |PlaneHeight               |
     |       |Height  |                          |
     |       |        +--------------------------+
     |       |        +--------------------------+ <- XaxisSubMargin
     |       |        |       PlaneNum = 0       |
     |       |        |PlaneHeight               |
     |       |        |                          |
     |       +        +------- PlaneWidth -------+
     |       |Xaxis
     |       |Margin
  0% +----------------+                          +---------+
        YaxisMargin                              RightMargin
%}
classdef MultiPlot < handle

% �萔�v���p�e�B
properties(Constant)

end

% �N���X�v���p�e�B
properties(SetAccess = protected)
	FigHandle;				% figure�n���h��
	FigLeft    = 0;			% [px] figure�̍����W
	FigBottom  = 0;			% [px] figure�̉����W
	FigWidth   = 0;			% [px] figure�̕�
	FigHeight  = 0;			% [px] figure�̍���
	NumPlanes  = 0;			% �v���b�g���ʂ̒i��
	CurrPlane  = 0;			% ���ݑI�𒆂̃v���b�g����
	PlaneHandle;			% axes�n���h��
	PlaneLeft  = 0;			% [px] �v���b�g���ʂ̍��ʒu
	PlaneWidth = 0;			% [px] �v���b�g���ʂ̉���
	PlaneHeight = 0;		% [px] �v���b�g���ʂ̍���
	AllPlanesHeight = 0;	% [px] ���ׂẴv���b�g���ʂ̍���
	YaxisMargin = 130;		% [px] Y���p�̗]��
	XaxisSubMargin = 3;		% [px] ���Ԓi��X���p�̗]��
	XaxisMargin = 70;		% [px] �ŉ��i��X���p�̗]��
	TopMargin   = 10;		% [px] �㑤�̗]��
	RightMargin = 20;		% [px] �E���̗]��
	FontSize   = 18;		% [pt] �����T�C�Y
	FontName   = 'Times New Roman';	% �t�H���g��
	AxisColor  = [0, 0, 0];			% ���̐F
	BackColor  = [1, 1, 1];			% �w�i�F
	ColorBlack = [      0,       0,       0];	% �v���b�g���̐F(���n)
	ColorGray  = [127/255, 127/255, 127/255];	% �v���b�g���̐F(�D�n)
	ColorRed   = [225/255,   0/255,   0/255];	% �v���b�g���̐F(�Ԍn)
	ColorGreen = [  0/255, 255/255, 127/255];	% �v���b�g���̐F(�Όn)
	ColorBlue  = [  0/255, 127/255, 255/255];	% �v���b�g���̐F(�n)
	WidthNarrow = 1.2;		% �v���b�g���̑���(��)
	WidthNormal = 2.0;		% �v���b�g���̑���(����)
	WidthWide   = 2.8;		% �v���b�g���̑���(��)
	LineColors;				% �v���b�g���̐F�p�A�z�z��
	LineWidths;				% �v���b�g���̑����p�A�z�z��
	LegendEdgeColor = [1, 1, 1];			% �}��̘g�F
	LegendBackColor = [0.95, 0.95, 0.95];	% �}��̔w�i�F
	XaxisName = 'Time [s]';	% ���Ԏ��̃��x����
	Xmin = 0;				% X���͈͂̍ŏ��l
	Xstp = 0;				% X���̃O���b�h�Ԋu
	Xmax = 0;				% X���͈͂̍ő�l
end

% �ÓI���\�b�h
methods(Static)
	% CSV�t�@�C����ǂݍ��ފ֐�
	function varargout = LoadCsvFile(FileName)
		% CSV�t�@�C���ǂݍ���
		FileFormat = FileName(length(FileName)-3:length(FileName));	% �g���q�ǂݍ���
		if strcmp(FileFormat, '.csv')
			RawData = readmatrix(FileName,'FileType','text');		% CSV�t�@�C���̏ꍇ
		else
			assert(false, 'MultiPlot: ���m�̃t�@�C���`��');			% ���m�̃t�@�C���̏ꍇ
		end
		
		% �o�͕ϐ��̐��̃`�F�b�N
		[h, w] = size(RawData);
		if w < nargout
			assert(false, 'MultiPlot: CSV�t�@�C���̗�̐� < �o�͕ϐ��̐�');
		end
		
		% �o�͕ϐ��Ƃ��ĕԂ�
		for i = 1:nargout
			varargout{i} = RawData(:, i);
		end
	end
end

% ���\�b�h
methods
	
	% �R���X�g���N�^
	function this = MultiPlot(Handle)
		this.FigHandle = Handle;
		clf(this.FigHandle);
		set(this.FigHandle, 'PaperPositionMode','auto');
		set(this.FigHandle, 'Position',[this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
		set(this.FigHandle, 'Color',this.BackColor);
		
		% �v���b�g���F�̏�����
		ColorDefs = {this.ColorBlack, this.ColorGray, this.ColorRed, this.ColorGreen, this.ColorBlue};
		ColorNames = ["Black", "Gray", "Red", "Green", "Blue"];
		this.LineColors = dictionary(ColorNames, ColorDefs);	% ���A�z�z��(R2022b�ȍ~���K�v)
		
		% �v���b�g�����̏�����
		WidthDefs = {this.WidthNarrow, this.WidthNormal, this.WidthWide};
		WidthSizes = ["Thin", "Normal", "Heavy"];
		this.LineWidths = dictionary(WidthSizes, WidthDefs);	% ���A�z�z��(R2022b�ȍ~���K�v)
	end
	
	% �O���t�ʒu�̐ݒ�֐�
	function FigurePosition(this, Left, Bottom)
		this.FigLeft = Left;
		this.FigBottom = Bottom;
		set(this.FigHandle, 'Position', [this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
	end
	
	% �O���t�T�C�Y�̐ݒ�֐�
	function FigureSize(this, Width, Height)
		this.FigWidth = Width;
		this.FigHeight = Height;
		set(this.FigHandle, 'Position', [this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
	end
	
	% �]���̐ݒ�֐�
	function FigureMargin(this, Left, Bottom, Right)
		this.YaxisMargin = Left;	% [px] Y���p�̗]��
		this.XaxisMargin = Bottom;	% [px] �ŉ��i��X���p�̗]��
		this.RightMargin = Right;	% [px] �E���̗]��
	end
	
	% �t�H���g�̐ݒ�֐�
	function Font(this, Name, Size)
		this.FontName = Name;
		this.FontSize = Size;
	end
	
	% �v���b�g���ʂ̒i���̐ݒ�֐�
	function NumOfPlanes(this, NumPlanes)
		% �T�C�Y��`�p�����[�^�̌v�Z
		this.NumPlanes = NumPlanes;
		this.AllPlanesHeight = this.FigHeight - this.XaxisMargin - this.TopMargin + this.XaxisSubMargin;% [px] ���ׂẴv���b�g���ʂ̍���
		this.PlaneHeight = this.AllPlanesHeight/this.NumPlanes - this.XaxisSubMargin;					% [px] �v���b�g���ʂ̍���
		this.PlaneWidth = this.FigWidth - this.YaxisMargin - this.RightMargin;							% [px] �v���b�g���ʂ̕�
		this.PlaneLeft = this.YaxisMargin;																% [px] �v���b�g���ʂ̍��ʒu
		
		% �i�����̃v���b�g���ʂ𐶐�
		for i = 0:1:(NumPlanes - 1)
			this.PlaneHandle(i + 1) = this.GeneratePlane(i);
		end
	end
	
	function FontSettings(this, Ax)
		% X�����x���ƃt�H���g�̐ݒ�
		if this.CurrPlane ~= 0
			set(Ax,'XTickLabel',{});	% �ŉ��i�ȊO��X�����x��������
		else
			xlabel(this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% �ŉ��i�̂�X�����x����t����
		end
		
		% �t�H���g�ƐF�̐ݒ�
		set(Ax, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'Color',this.BackColor, 'XColor',this.AxisColor, 'YColor',this.AxisColor);
	end
	
	% �v���b�g���ʂ𐶐�����֐�
	function Ax = GeneratePlane(this, PlaneNum)
		% ���W�v�Z
		w = 1/this.FigWidth*this.PlaneWidth;	% �v���b�g���ʂ̕������Z [px] �� [%]
		h = 1/this.FigHeight*this.PlaneHeight;	% �v���b�g���ʂ̍��������Z [px] �� [%]
		l = 1/this.FigWidth*this.PlaneLeft;		% �v���b�g���ʂ̍��ʒu�����Z [px] �� [%]
		b = 1/this.FigHeight*(this.AllPlanesHeight/this.NumPlanes*PlaneNum + this.XaxisMargin);		% �v���b�g���ʂ̉��ʒu�����Z [px] �� [%]
		Pos = [l, b, w, h];						% [��%, ��%, ��%, ����%] 
		Ax = axes('Position', Pos);				% �v���b�g���ʐ���
		plot(Ax, 0, 0);							% �������pplot() ��v���b�g�����Ă���
		
		% X�����x���ƃt�H���g�̐ݒ�
		if PlaneNum ~= 0
			set(Ax,'XTickLabel',{});	% �ŉ��i�ȊO��X�����x��������
		else
			xlabel(this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% �ŉ��i�̂�X�����x����t����
		end
		
		% �t�H���g�ƐF�̐ݒ�
		ylabel(Ax, 'y', 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'Color',this.BackColor, 'XColor',this.AxisColor, 'YColor',this.AxisColor);
	end
	
	% �v���b�g���ʂ�I������֐�(���i�ڂɃv���b�g���邩�̐ݒ�)
	function SelectPlane(this, PlaneNum)
		this.CurrPlane = this.NumPlanes - PlaneNum;
		axes(this.PlaneHandle( this.NumPlanes - PlaneNum + 1 ))
	end
	
	% X���͈͂̐ݒ�֐�
	function XaxisRange(this, Xmin, Xstp, Xmax)
		this.Xmin = Xmin;
		this.Xstp = Xstp;
		this.Xmax = Xmax;
		
		% �i�����̃v���b�g���ʂɑ΂���X���͈͂�ݒ�
		for i = 1:1:this.NumPlanes
			xlim(this.PlaneHandle(i), [Xmin, Xmax]);
			set(this.PlaneHandle(i),'XTickMode','manual','XTick', Xmin:Xstp:Xmax);
		end
	end
	
	% X�����x�����̐ݒ�֐�
	function XaxisLabel(this, Xname)
		this.XaxisName = Xname;
		xlabel(this.PlaneHandle(1), this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% �ŉ��i��X�����x��
	end
	
	% X/Y���O���b�h�������ݒ肷��֐�
	function AutoGrid(this, Ax)
		Ymax = max(Ax.YTick);
		Ymin = min(Ax.YTick);
		Ystp = (Ymax - Ymin)/(length(Ax.YTick) - 1);
		Xmax = max(Ax.XTick);
		Xmin = min(Ax.XTick);
		Xstp = (Xmax - Xmin)/(length(Ax.XTick) - 1);
		grid on;
		axis([Xmin, Xmax, Ymin - Ystp/2, Ymax + Ystp/2]);
		set(Ax,'XTickMode','manual','YTickMode','manual','XTick', Xmin:Xstp:Xmax,'YTick', Ymin:Ystp:Ymax);
		Ax.YAxis.Exponent = 0;
	end
	
	% Y���O���b�h���蓮�ݒ肷��֐�
	function ManualGrid(this, Ymin, Ystp, Ymax)
		grid on;
		ylim([Ymin - Ystp/2, Ymax + Ystp/2]);
		set(gca,'YTickMode','manual','YTick', Ymin:Ystp:Ymax);
	end
	
	% �c�����x���̐ݒ�֐�
	function Label(this, Name)
		ylabel(gca, Name, 'FontSize',this.FontSize, 'FontName',this.FontName);
	end
	
	% �}��̐ݒ�֐�
	function Legend(this, Name, Position, Direction)
		nll = {char};		% �󕶎�
		Name = [nll, Name];	% �������pplot()�̕��̖}��͔�\��
		h = legend(Name, 'Location',Position, 'Orientation',Direction);
		set(h, 'TextColor',this.AxisColor);
		set(h, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(h, 'EdgeColor', this.LegendEdgeColor, 'Color', this.LegendBackColor);
	end
	
	% �f�[�^���_�E���T���v�����O���ĊԈ����֐�
	function xout = DownSampling(this, xin, dsmpl)
		len = length(xin(:,1));
		xout = xin(1:dsmpl:len, :);
	end
	
	% �v���b�g�֐�
	function Plot(this, varargin)
		% ������ǂݍ���
		PlotFunc = varargin{1};	% �֐��n���h��
		Argin = varargin{2};	% �����Z���z��
		x = Argin{1};			% X���f�[�^
		y = Argin{2};			% Y���f�[�^
		Color = Argin{3};		% ���̐F
		Width = Argin{4};		% ���̑���
		if size(varargin) == [1, 3]
			% �v���b�g�_�̎�ނ��w�肳�ꂽ�ꍇ
			PlotType = varargin{3};
		else
			% �v���b�g�_�̎�ނ��w�肳��Ȃ������ꍇ
			PlotType = ' ';
		end
		
		% �����̐��œ����ύX
		if length(Argin) == 4
			% �_�E���T���v�����O���Ȃ��ꍇ
			DownSmplRate = 1;
		else
			% �_�E���T���v�����O����ꍇ
			DownSmplRate = Argin{5};
			x = this.DownSampling(x, DownSmplRate);
			y = this.DownSampling(y, DownSmplRate);
		end
		
		% ���ۂ̃v���b�g����
		hold on;
		h = PlotFunc(gca, x, y, PlotType);
		c = this.LineColors(Color);
		w = this.LineWidths(Width);
		set(h, 'LineWidth', w{1}, 'Color', c{1});
		hold off;
		this.AutoGrid(gca);
	end
	
	% �K�i�v���b�g�֐�
	function StairsPlot(this, varargin)
		this.Plot(@stairs, varargin);
	end
	
	% ���`�v���b�g�֐�
	function LinePlot(this, varargin)
		this.Plot(@plot, varargin);
	end
	
	% �_�v���b�g�֐�
	function DotPlot(this, varargin)
		this.Plot(@plot, varargin, '.');
	end
	
	% �~�v���b�g�֐�
	function CrossPlot(this, varargin)
		this.Plot(@plot, varargin, 'x');
	end
	
	% ���v���b�g�֐�
	function CirclePlot(this, varargin)
		this.Plot(@plot, varargin, 'o');
	end
	
	% �{�v���b�g�֐�
	function PlusPlot(this, varargin)
		this.Plot(@plot, varargin, '+');
	end
	
	% PNG�摜��EPS�t�@�C������
	function SavePNGandEPS(this, FileName)
		PngFileName = strcat(FileName(1:length(FileName)-4), '.png');	% �g���q��PNG�ɕϊ�
		EpsFileName = strcat(FileName(1:length(FileName)-4), '.eps');	% �g���q��EPS�ɕϊ�
		exportgraphics(this.FigHandle, PngFileName,'ContentType','image', 'Resolution',150);
		exportgraphics(this.FigHandle, EpsFileName,'ContentType','vector','BackgroundColor','none');
	end
	
	% PNG�摜��PDF�t�@�C������
	function SavePNGandPDF(this, FileName)
		PngFileName = strcat(FileName(1:length(FileName)-4), '.png');	% �g���q��PNG�ɕϊ�
		PdfFileName = strcat(FileName(1:length(FileName)-4), '.pdf');	% �g���q��PDF�ɕϊ�
		exportgraphics(this.FigHandle, PngFileName,'ContentType','image', 'Resolution',150);
		exportgraphics(this.FigHandle, PdfFileName,'ContentType','vector');
	end
end
end

