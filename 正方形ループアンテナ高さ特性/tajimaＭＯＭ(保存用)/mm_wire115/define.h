//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			�ϐ��A�萔�A�z��錾�w�b�_�[�@�@define.h
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//				�v���g�^�C�v�錾
//////////////////////////////////////////////////////////////////////////////////
//--������
void Initialization(void);						//�z��ƌW���̏�����
//--�z��
void MakeMatAll(void);							//�z��m��
void DelMatAll(void);							//�z��J��
//--�v�Z
void Calculation(void);							//�v�Z����
//--�C���s�[�_���X�s��֌W
void MakeZmn(void);								//�C���s�[�_���X�s��쐬
void MakeGaussPara(void);						//�K�E�X�ϕ��̃p�����[�^�쐬
complex grmn(int ,int ,int ,int);				//g��mn�̌v�Z (m ,n ,i ,j)
complex hij(int ,int ,int ,int ,int ,double);	//hji (<1 or 2 or 3>,m ,n ,i ,j ,xi)
//--�d�����z�֌W
void MakeCurrent(void);							//�d�����z
void MakePhase(void);							//�d���ʑ��v�Z
void Pibot(int);								//�s�{�b�g�̑I��(Zmn[i][i]��0�̔���)
void Gauss(void);								//�K�E�X�̏����@
complex Idx(int ,int , double);					//�d���v�Z (���C���[No 
												//		   ,���C���[���̃Z�O�����gNo
												//         ,�ʒu)
//--���ˊE�֌W
void Radiation(void);							//���ˊE�v�Z
complex FarFieldT(double ,double);				//�����E�Ɛ��� (�ƁC��)����
complex FarFieldF(double ,double);				//�����E�Ӑ��� (�ƁC��)����
//--���o�͊֌W
void ConfInput(void);							//����
void OutputCurrent(void);						//�d���o��
void OutputConf(void);							//�`��o��
void OutputRad(void);							//���ˊE�@���˃p�^�[���o��
void OutputChara(void);							//�����f�[�^�o��
void FileOpen(void);							//�t�@�C���J��
void FileClose(void);							//�t�@�C������
void StartTime(void);							//�^�C�}�[�X�^�[�g
void LapsedTime(void);							//�o�ߎ���
void ConfCheck(void);							//�`�󌟍�

//////////////////////////////////////////////////////////////////////////////////
//			�萔
//////////////////////////////////////////////////////////////////////////////////
//--��{�����萔
const double PI = 3.1415926535897932385;		//��
const double C  = 2.99792458*pow(10.0,8.0);			//����
const double e0 = 8.85418780*pow(10.0,-12.0);		//�^�󒆂̓�������0
const double u0 = 4.0*PI*pow(10.0,-7.0);			//�^�󒆂̗U�d����0
const double R  = 1.0;							//1.0 �Œ�(�����E�v�Z����)
//--���ȗ��p
complex J=Complex(0.0,1.0);                     //�����P�ʁ@J=��-1

//--�ϕ��p
const int GaussTenNor = 4 ;						//�K�E�X�ϕ����_  [ 4�Œ�]�i�ʏ�p�j
const int GaussTenSpe = 40;						//�K�E�X�ϕ����_��[40�Œ�]�i���ٓ_�p�j

//////////////////////////////////////////////////////////////////////////////////
//			�z��̖���
//////////////////////////////////////////////////////////////////////////////////
//--�`��z��
double *RX , *RY , *RZ;							//�e�Z�O�����g�̎n�_�@	RX[ ]  RY[ ]  RZ[ ]
double *SX , *SY , *SZ;							//�P�ʃx�N�g���̐����@	SX[ ]  SY[ ]  SZ[ ]
double *SEGL , *RA;								//�e�Z�O�����g�̒��� ���a	SEGL[ ]	RA[ ]�@		
int *SEGN , *FEDP;								//�e���C���[�̃Z�O�����g���@SEGN[ ]
complex *FEDV;									//�e���d�d��			FEDV[ ]
int *RAD_MODE;									//�e�Z�O�����g�̕��˃��[�h		RAD_MODE[ ]
int *LOADP;										//�e�C���s�[�_���X���׈ʒu
complex *LOADZ;									//�e�C���s�[�_���X���גl

//--�d���z��
complex **Zmn;									//�C���s�[�_���X�s��	Zmn[ ][ ]
complex *Im;									//�d�����z				Im[ ]
double *PhaseIm;								//�d���̈ʑ�			PhaseIm[ ]

//--���̓C���s�[�_���X
complex *Zin;									//�e���d�_�̓��̓C���s�[�_���X
double *VSWR_50 , *VSWR_75;						//�e���d�_��50����75���ɑ΂���VSWR

//--���ˊE�z��
complex *TRAD ,*FRAD;							//E�Ƃ�E�ӂ̕��ˊE TRAD[ ] FRAD[ ]
complex *RRAD ,*LRAD;							//ER ��EL �̕��ˊE RRAD[ ] LRAD[ ]
double *TPHASE ,*FPHASE;						//E�Ƃ�E�ӂ̈ʑ� TPHASE[ ] FPHASE[ ]
double *RPHASE ,*LPHASE;						//ER ��EL �̈ʑ� RPHASE[ ] LPHASE[ ]
double *TGAI ,*FGAI;							//E�Ƃ�E�ӂ̗��� TGAI[ ] FGAI[ ]
double *TFGAI;									//E�Ƃ�E�ӂ̍��v�̗��� TFGAI[ ]
double *RLGAI;									//ER ��EL �̍��v�̗��� RLGAI[ ]
double *RGAI ,*LGAI;							//ER ��EL �̗��� RGAI[ ] LGAI[ ]
double *TPAT ,*FPAT;							//E�Ƃ�E�ӂ̃p�^�[�� TPAT[ ] FPAT[ ]
double *RPAT ,*LPAT;							//ER ��EL �̃p�^�[�� RPAT[ ] LPAT[ ]
double *AR;										//���� AR[ ]

//--�K�E�X�ϕ��̃p�����[�^
double *GaussWeightNor , *GaussBuntenNor;		//�d�� ���_
double *GaussWeightSpe , *GaussBuntenSpe;		//�d�� ���_ ���ٓ_�p

//////////////////////////////////////////////////////////////////////////////////
//			�t�@�C���|�C���^
//////////////////////////////////////////////////////////////////////////////////
FILE *fp_conf;									//�`��o�͗p
FILE *fp_curr;									//�d�����z�o�͗p
FILE *fp_radf;									//���ˊE�o�͗p
FILE *fp_radp;									//���˃p�^�[���o�͗p
FILE *fp_char1;									//�����o�͗p
FILE *fp_char2;									//�����o�͗p
FILE *fp_free;									//�C�Ӄf�[�^�o�͗p

//////////////////////////////////////////////////////////////////////////////////
//			�ϐ�
//////////////////////////////////////////////////////////////////////////////////
//--�`�� ���d
double FREQ0;									//�݌v���g��f0
double RAMDA0;									//�݌v���g���̎��R��Ԕg��
double USEF;									//���d���g��
int NWIR;										//�S���C���[��
int NSEG;										//�S�Z�O�����g��
int NFED;										//�S���d�_��
int NSEG0;										//�S�d���v�Z�_��
int NLOAD;										//��R�����א�

//--�v�Z�p
double k0;										//k0
complex uair;									//u(air)
double Beta;									//��

//--���ˊE�v�Z
static int    AXMODE   ;						//(0=�ӌŒ�, 1=�ƌŒ�)	
double DegDelta ;								//���ݕ�("0.0"�֎~)
double DegStart ;								//�����p
double DegWidth ;								//�͈�
double FixAngle ;								//�Œ莲�p�x

//--���ԑ���
int		la_min,la_sec;							//�� �b
time_t	start, finish;							//�J�n�A�I��

//--�v�Z�񂵗p�p�����[�^
double PARA1,PARA2,PARA3;	//�o�͂Ɏg������