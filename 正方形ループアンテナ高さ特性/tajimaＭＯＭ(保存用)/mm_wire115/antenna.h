//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			�A���e�i�`��w�b�_�[�@�@antenna.h
//				�_�C�|�[���A���e�i
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//				�v�Z
//////////////////////////////////////////////////////////////////////////////////
void Calculation(void)
{
	//========================================================
	//			���[�J���ϐ�
	//========================================================
	//--�J�E���^
	int i,j,l,o,p,f,n_seg_free;

	//--�v�Z�Ŏg�p����ϐ�
	double rrx , rry , rrz;				//�Z�O�����g�n�_
	double drx , dry , drz;				//�Z�O�����g�I�_
	double ssx , ssy , ssz;				//�Z�O�����g�̒P�ʃx�N�g��
	double dcl;							//�Z�O�����g����

	//--�G�������g���֌W�i���C���[�{���j
	int n_array;						//�_�C�|�[���̖{��

	//--�e�����̕�����
	int nbun_pole1,nbun_pole2,nbun_pole3,n_pole;//�_�C�|�[���̕�����
	
	//--�e���̃Z�O�����g��
	double r_pole1,r_pole2,r_pole3,r_pole;//�_�C�|�[���̃Z�O�����g��

	//--�݌v�l�֌W
	double L, Lhalf, DeruL, h;//�_�C�|�[���̒���
	double wire_radius;					//���C���[���a

	//double optim[5][50][2];

	////s=0�̂Ƃ��]���A���e�i
	//optim[0][0][0] = 0.254;/*0.1_L*/ optim[0][0][1] = 0.105;/*0.1_dL*/
	//optim[0][1][0] = 0.253;/*0.11_L*/ optim[0][1][1] = 0.112;/*0.11_dL*/
	//optim[0][2][0] = 0.251;/*0.12_L*/ optim[0][2][1] = 0.120;/*0.12_dL*/
	//optim[0][3][0] = 0.249;/*0.13_L*/ optim[0][3][1] = 0.126;/*0.13_dL*/
	//optim[0][4][0] = 0.247;/*0.14_L*/ optim[0][4][1] = 0.132;/*0.14_dL*/
	//optim[0][5][0] = 0.246;/*0.15_L*/ optim[0][5][1] = 0.136;/*0.15_dL*/
	//optim[0][6][0] = 0.244;/*0.16_L*/ optim[0][6][1] = 0.141;/*0.16_dL*/
	//optim[0][7][0] = 0.242;/*0.17_L*/ optim[0][7][1] = 0.145;/*0.17_dL*/
	//optim[0][8][0] = 0.242;/*0.18_L*/ optim[0][8][1] = 0.146;/*0.18_dL*/
	//optim[0][9][0] = 0.240;/*0.19_L*/ optim[0][9][1] = 0.150;/*0.19_dL*/
	//optim[0][10][0] = 0.238;/*0.2_L*/ optim[0][10][1] = 0.153;/*0.20_dL*/
	//optim[0][11][0] = 0.237;/*0.21_L*/ optim[0][11][1] = 0.155;/*0.21_dL*/
	//optim[0][12][0] = 0.237;/*0.22_L*/ optim[0][12][1] = 0.155;/*0.22_dL*/
	//optim[0][13][0] = 0.236;/*0.23_L*/ optim[0][13][1] = 0.157;/*0.23_dL*/
	//optim[0][14][0] = 0.235;/*0.24_L*/ optim[0][14][1] = 0.159;/*0.24_dL*/
	//optim[0][15][0] = 0.235;/*0.25_L*/ optim[0][15][1] = 0.159;/*0.25_dL*/
	

	//========================================================
	//			���ˊE�v�Z�ݒ�		�P�ʁF(deg)
	//========================================================
	//--�v�Z���[�h
	AXMODE = 0;			//(0=�ӌŒ�, 1=�ƌŒ�)

	//--�ώ�
	DegDelta =    1.0;	//���ݕ�("0.0"�֎~)
	DegStart =  0.0;	//�����p
	DegWidth =  360.0;	//�͈�

	//--�Œ莲
	FixAngle =  0.0;	//�Œ�p�x

	//========================================================
	//			�݌v���g��		FREQ0	(f0)	�P�ʁF(Hz)
	//			���d���g��		USEF	(f)		�P�ʁF(Hz)
	//========================================================
	//--�݌v���g���ݒ�
	//for(k=0.17;k<=0.25;k+=0.01){

		//for(m=0.17;m<0.25;m+=0.01){//�ۓ��f�q
			
			//for(n=-0;n<1;n++){//���g��

	for(l=10;l<20;l++){//L

		for(o=0;o<1;o++){//����

			for(p=0;p<1;p++){//dL

				for(f=0;f<1;f++){//���g��





	//--�݌v���g���ݒ�@
	FREQ0 = 1.0 * pow(10.0,9.0);

	//--�d�������ݒ�@
	USEF  = FREQ0 * (1.000 + 0.001*f);				//f���Ƃ�ɂ͌W����ύX

	//--�݌v���g���̎��R��Ԕg��
	RAMDA0= C/FREQ0;

	int s = 0;

	//========================================================
	//			�`��l��`
	//			������RAMDA0����@�W���Œl��`(������)
	//========================================================
	//--����
	L = (0.10 + 0.01 * l) * RAMDA0;
	Lhalf = L/2;
	DeruL = (0.05 + 0.01 * p) * RAMDA0;
	h = (0.10 + 0.01*o) * RAMDA0;
	
	r_pole1 = Lhalf;
	r_pole2 = DeruL;   
	r_pole3 = h;


	//--���C���[���a
	wire_radius = RAMDA0*0.001;

	//========================================================
	//			�S���C���[��	NWIR
	//			�Z�O�����g��	NSEG
	//			���d�_��		NFED
	//			�S�d���v�Z�_��	NSEG0
	//			�y���א�		NLOAD		
	//========================================================
	//--�S���C���[���ݒ�
	n_array = 19;
	NWIR = n_array;

	//--����������
	nbun_pole1 = int(r_pole1/(0.025*RAMDA0));//a
	nbun_pole2 = int(r_pole2/(0.025*RAMDA0));//b
	nbun_pole3 = int(r_pole3/(0.025*RAMDA0));//c

	//--�S�Z�O�����g���ݒ�
	n_pole = 2 * (nbun_pole1*8+8 + nbun_pole2+1) + nbun_pole3+2;					
	NSEG = n_pole;						//�S�Z�O�����g��

	//printf("%d\n",NSEG);

	//--�S�d���v�Z�_��  [NSEG0=NSEG-NWIR]	�������������폜�@�s��
	NSEG0 = NSEG - NWIR;

	//--���d�_��
	NFED = 1;

	//--�C���s�[�_���X���א�
	NLOAD = 0;

	//========================================================
	//			�z��m�ۂƏ�����
	//			[���ˊE�v�Z�ݒ�,NWIR,NSEG,NFED,NLOAD]
	//				�ݒ��Ɋ֐��Ăяo��
	//			�������������폜�@�s��
	//========================================================
	MakeMatAll();								//�z��m��(���R�����g�A�E�g�֎~��)
	Initialization();							//�z��ƌv�Z�p�ϐ��̏�����



	//========================================================
	//			�o�̓t�@�C���p�̃p�����[�^
	//			PARA1 �� PARA2 �ɒl����͂���D
	//========================================================
	//PARA1 = r_pole/RAMDA0;						//�o�͗p�p�����[�^�@
	//PARA2 = wire_radius;						//�o�͗p�p�����[�^�A
	PARA1 = L/RAMDA0;									//�o�͗p�p�����[�^�B
	PARA2 = DeruL/RAMDA0;								//�o�͗p�p�����[�^�C
	PARA3 = h/RAMDA0;									//�o�͗p�p�����[�^�D

	//========================================================
	//			�݌v�l����		�P�ʁF(m)
	//			�Z�O�����g�n�_	RX[ ] , RY[ ] , RZ[ ] 
	//			�P�ʃx�N�g��	SX[ ] , SY[ ] , SZ[ ]
	//			�Z�O�����g��	SEGL[ ]
	//========================================================
i=0;
//----------------------------------------------------------------------------------------------------


//�����̉���------------------------------------------------------------------------------------------
for(j = 0; j < nbun_pole1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (0-j) * (r_pole1/nbun_pole1);			//�n�_�ʒu�v�Zx
		rry = -r_pole1;											//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = (-1-j) * (r_pole1/nbun_pole1);			//�I�_�ʒu�v�Zx
		dry = -r_pole1;											//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		//printf("%d\n",SEGL[0]);
		i++;
	}

//�I�[�o�[���b�v------------------------------------------------------------------------------------------
for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = (j-nbun_pole1) * (r_pole1/nbun_pole1);			//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = (j-nbun_pole1+1) * (r_pole1/nbun_pole1);			//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}

//�����̏c��------------------------------------------------------------------------------------------
for(j = 0; j < nbun_pole1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = (j-nbun_pole1) * (r_pole1/nbun_pole1);			//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = (j-nbun_pole1+1) * (r_pole1/nbun_pole1);			//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}


//�I�[�o�[���b�v------------------------------------------------------------------------------------------
for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = (j+0) * (r_pole1/nbun_pole1);						//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = (j+1) * (r_pole1/nbun_pole1);						//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}

//���̏㔼��------------------------------------------------------------------------------------------
for(j = 0; j < nbun_pole1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = (j+0) * (r_pole1/nbun_pole1);						//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = (j+1) * (r_pole1/nbun_pole1);						//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}


//�I�[�o�[���b�v------------------------------------------------------------------------------------------
for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (j-nbun_pole1) * (r_pole1/nbun_pole1);			//�n�_�ʒu�v�Zx
		rry = r_pole1;											//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz
		drx = (j-nbun_pole1+1) * (r_pole1/nbun_pole1);			//�I�_�ʒu�v�Zx
		dry = r_pole1;											//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}

//����̉���------------------------------------------------------------------------------------------
for(j = 0; j < nbun_pole1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (j-nbun_pole1) * (r_pole1/nbun_pole1);			//�n�_�ʒu�v�Zx
		rry = r_pole1;											//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = (j-nbun_pole1+1) * (r_pole1/nbun_pole1);			//�I�_�ʒu�v�Zx
		dry = r_pole1;											//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}


//�I�[�o�[���b�v------------------------------------------------------------------------------------------
for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (j+0) * (r_pole1/nbun_pole1);						//�n�_�ʒu�v�Zx
		rry = r_pole1;											//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz
		drx = (j+1) * (r_pole1/nbun_pole1);						//�I�_�ʒu�v�Zx
		dry = r_pole1;											//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}


//�E���̔��˔�------------------------------------------------------------------------------------------

	for(j = 0; j<nbun_pole1*4+1+1+1+1; ++j){
		//--�v�Z���ʑ��
		RX[i] = -RX[j];											//�n�_�ʒux
		RY[i] = -RY[j];											//�n�_�ʒuy
		RZ[i] = 0;												//�n�_�ʒuz
		SX[i] = -SX[j];											//�P�ʃx�N�g��x
		SY[i] = -SY[j];											//�P�ʃx�N�g��y
		SZ[i] = 0;												//�P�ʃx�N�g��z
		SEGL[i] = SEGL[j];										//�Z�O�����g��
		i++;
	}
//--------------------------------------------------------------------------------------------------


//�I�[�o�[���b�v�ۓ��f�q------------------------------------------------------------------------------------------
for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = 0;												//�n�_�ʒu�v�Zx
		rry = (j-nbun_pole1) * (r_pole2/nbun_pole2);			//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = 0;												//�I�_�ʒu�v�Zx
		dry = (j-nbun_pole1+1) * (r_pole2/nbun_pole2);			//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}

double temp=0;
temp = dry;

//�ۓ��f�q------------------------------------------------------------------------------------------
for(j = 0; j < nbun_pole2; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = 0;												//�n�_�ʒu�v�Zx
		rry = (j-nbun_pole1) * (r_pole2/nbun_pole2);			//�n�_�ʒu�v�Zy
		rrz = 0;												//�n�_�ʒu�v�Zz																								
		drx = 0;												//�I�_�ʒu�v�Zx
		dry = (j-nbun_pole1+1) * (r_pole2/nbun_pole2);			//�I�_�ʒu�v�Zy
		drz = 0;												//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
		i++;
	}

//�S�Ă̔��˔�------------------------------------------------------------------------------------------

	for(j = 0; j<2*(nbun_pole1*4+1+1+1+1)+nbun_pole2+1; ++j){
		//--�v�Z���ʑ��
		RX[i] = RX[j];											//�n�_�ʒux
		RY[i] = RY[j];											//�n�_�ʒuy
		RZ[i] = -r_pole3;									//�n�_�ʒuz
		SX[i] = SX[j];											//�P�ʃx�N�g��x
		SY[i] = SY[j];											//�P�ʃx�N�g��y
		SZ[i] = -r_pole3;									//�P�ʃx�N�g��z
		SEGL[i] = SEGL[j];										//�Z�O�����g��
		i++;
	}
//--------------------------------------------------------------------------------------------------


//------------------------------------------------�G�I�[�o�[���b�v:���d�_(-z)------------------------------------------------------------
	for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = -r_pole1;											//�n�_�ʒu�v�Zy ���a
		rrz = -r_pole1;		//�n�_�ʒu�v�Zz
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = -r_pole1;												//�I�_�ʒu�v�Zy ���a
		drz = (j-1) * (r_pole3/nbun_pole3);		//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;
		i++;
	}

//------------------------------------------------�C���d�_------------------------------------------------------------

	for(j = 0; j < nbun_pole3; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = -r_pole1;											//�n�_�ʒu�v�Zy 
		rrz = (j-nbun_pole3) * (r_pole3/nbun_pole3);		//�n�_�ʒu�v�Zz
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = -r_pole1;												//�I�_�ʒu�v�Zy ���a
		drz = (j-nbun_pole3+1) * (r_pole3/nbun_pole3);		//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;
		i++;
	}
	

//------------------------------------------------�G�I�[�o�[���b�v:���d�_(+z)------------------------------------------------------------
	for(j = 0; j < 1; ++j){
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -r_pole1;											//�n�_�ʒu�v�Zx
		rry = -r_pole1;											//�n�_�ʒu�v�Zy ���a
		rrz = -(j+0) * (r_pole3/nbun_pole3);		//�n�_�ʒu�v�Zz
		drx = -r_pole1;											//�I�_�ʒu�v�Zx
		dry = -r_pole1;												//�I�_�ʒu�v�Zy ���a
		drz = -(j+1) * (r_pole3/nbun_pole3);		//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;
		i++;
	}

	//========================================================
	//			�e�Z�O�����g�̃��C���[���a����	�P�ʁF(m)
	//			RA[ ]		���C���[���a
	//========================================================
	for(i = 0; i < NSEG; ++i){
		RA[i] = wire_radius;
	}

	//========================================================
	//			�e�G�������g(���C���[)���̃Z�O�����g��
	//			SEGN[ ]
	//========================================================
//*	
	SEGN[0] = nbun_pole1 + 1;
    SEGN[1] = nbun_pole1 + 1;
	SEGN[2] = nbun_pole1 + 1;
	SEGN[3] = nbun_pole1 + 1;
	SEGN[4] = nbun_pole1 + 1;
    SEGN[5] = nbun_pole1 + 1;
	SEGN[6] = nbun_pole1 + 1;
	SEGN[7] = nbun_pole1 + 1;
	SEGN[8] = nbun_pole2 + 1;
	SEGN[9] = nbun_pole1 + 1;
    SEGN[10] = nbun_pole1 + 1;
	SEGN[11] = nbun_pole1 + 1;
	SEGN[12] = nbun_pole1 + 1;
	SEGN[13] = nbun_pole1 + 1;
    SEGN[14] = nbun_pole1 + 1;
	SEGN[15] = nbun_pole1 + 1;
	SEGN[16] = nbun_pole1 + 1;
	SEGN[17] = nbun_pole2 + 1;
	SEGN[18] = nbun_pole3 + 2;

	//*/
	//========================================================
	//			���d�֌W�ݒ�
	//			���d�ʒu		FEDP[ ]
	//			���d�d��		FEDV[ ]
	//
	//			���ʒu  FEDP[ ]
	//					[���d�������ʒu]-[���d���������C���[(���Ԗ�)]
	//			���d��  FEDV[ ]
	//				�@��P���d�̓d��
	//					�s���t���d��p����ꍇ   -2�~(1+j0)[V]
	//					���t���d��p����ꍇ     -1�~(1+j0)[V]
	//				�A���˔t���_�C�|�[����(���˔t�����t���d�A���e�i)
	//					��P���d(���˔\��)   -1�~(1+j0)[V]
	//					��Q���d(���˔���)    1�~(1+j0)[V]
	//						(���˔̕\�Ɨ��ŋ��d�ʑ����t�ɂ���)
	//				�B��P���d���C���[�W�@�̏ꍇ
	//					���̋��d�_���C���[�W�@��K�p����D
	//				�C��P���d���C���[�W�@��p���Ȃ��ꍇ
	//					���̋��d�_���C���[�W�@��p���Ȃ��D
	//				�D�������d�̓d���E�ʑ�������
	//					��P���d�͏����@�A��K�p
	//					��P���d�ȊO�Ő�����s���D
	//========================================================
	FEDP[0] = SEGN[0] - 1 - 1;		// ���d�ʒu
	FEDV[0] = -2*Complex(1.0 , 0.0);	// 1.0+j0.0[V]

//n_seg_free=0;
//for(j=0;j<n_array;j++){
//printf("SEGN[%d]=%d\n",j,SEGN[j]);
//n_seg_free+=SEGN[j];
//}
//printf("n_array=%d\n",n_array);
//printf("n_seg_free=%d\n",n_seg_free);
//printf("n_pole=%d\n",n_pole);


	//========================================================
	//			�C���s�[�_���X���׊֌W�ݒ�
	//			���׈ʒu		LOADP[ ]
	//			���גl�@		LOADZ[ ]
	//
	//			���ʒu    LOADP[ ]
	//					[���ׂ������ʒu]-[���ׂ��������C���[(���Ԗ�)]
	//			�����גl�@LOADZ[ ]  (50����ڑ�����ꍇ)
	//				�@���˔ɒ�R��t����ꍇ		 -2�~(50+j0)[��]
	//				�A���C���̓r���ɒ�R��t����ꍇ -1�~(50+j0)[��]
	//				�B�ڑ��Ώۂɂ�炸�C�l�́u�}�C�i�X(-)�v�ɂȂ�
	//								(��R�l�̓x�N�g���ł͖�������)
	//========================================================


	//========================================================
	//			�e�Z�O�����g����̕��ˏ�Ԃ�ݒ�(�������ːݒ�)
	//			RAD_MODE[ ] = 1 or 0
	//				(0 : ���˃J�b�g   1 : �ʏ�̕���)
	//			�������őS��RAD_MODE[ ]=1 �ɐݒ�ς�
	//========================================================


	//========================================================
	//			�d���E���ˊE�v�Z
	//========================================================	
	ConfCheck();				//�`�󌟍�
	MakeZmn();					//Zmn�쐬
	MakeCurrent();				//�d���v�Z
	MakePhase();				//�d���ʑ��v�Z
	Radiation();				//���ˊE�v�Z
	LapsedTime();				//�o�ߎ��Ԏ擾

	//========================================================
	//			�f�[�^�o��
	//========================================================	
	OutputConf();				//�`��o��
	OutputCurrent();			//�d���o��
	OutputRad();				//���ˊE�o��
	OutputChara();				//�����f�[�^�o��

	//========================================================
	//			�C�Ӄf�[�^�o��	free.dat
	//				�g�p�t�@�C���|�C���^ [fp_free]
	//========================================================	
	fprintf(fp_free," PARA1 = %9.5f    PARA2 = %9.5f   PARA3 = %9.5f\n",PARA1,PARA2,PARA3);	//�p�����^�[
	fprintf(fp_free," Zin = %9.2f %9.2f\n",Zin[0].re,Zin[0].im);		//�C���s�[�_���X
	fprintf(fp_free," �v�Z���� = %d:%d\n",la_min,la_sec);				//�o�ߎ��ԕ\��
	fprintf(fp_free," \n\n\n");											//���s

	//========================================================
	//			�z��̊J��
	//========================================================	
	DelMatAll();				//�z��J��(���R�����g�A�E�g�֎~��)
}
}
}
}

}


//============================================================================
//					�ϐ��Ɣz��ꗗ
//============================================================================
//
//--�`�� ���d
//		double	FREQ0					�݌v���g��f0
//		double	RAMDA0					�݌v���g���̎��R��Ԕg��
//		double	USEF					���d���g��
//		int		NWIR					�S���C���[��
//		int		NSEG					�S�Z�O�����g��
//		int		NFED					�S���d�_��
//		int		NSEG0					�S�d���v�Z�_��
//		int		NLOAD					��R�����א�
//--���ˊE�v�Z
//		int		AXMODE   				(0=�ӌŒ�, 1=�ƌŒ�)	
//		double	DegDelta 				���ݕ�
//		double	DegStart 				�����p
//		double	DegWidth 				�͈�
//		double	FixAngle 				�Œ莲�p�x
//--���ԑ���
//		int		la_min	la_sec			�� �b
//--�v�Z�񂵗p�p�����[�^
//		double	PARA1,PARA2				�p�����[�^�@�ƇA
//--�`��z��
//		double	RX[ ]  RY[ ]  RZ[ ]		�e�Z�O�����g�̎n�_
//		double�@SX[ ]  SY[ ]  SZ[ ]		�P�ʃx�N�g���̐���
//		double	SEGL[ ]	RA[ ]			�e�Z�O�����g�̒��� ���a�@		
//		int�@	SEGN[ ]					�e���C���[�̃Z�O�����g��
//		int		FEDP[ ]					�e���d�ʒu
//		complex	FEDV[ ]					�e���d�d��
//		int		RAD_MODE[ ]				�e�Z�O�����g�̕��˃��[�h
//		int		LOADP[ ]				�e�C���s�[�_���X���׈ʒu
//		complex LOADZ[ ]				�e�C���s�[�_���X���גl
//--�d���z��
//		complex Zmn[ ][ ]				�C���s�[�_���X�s��	
//		complex Im[ ]					�d�����z
//		double  PhaseIm[ ]				�d���̈ʑ�
//--���̓C���s�[�_���X
//		complex Zin[ ]					�e���d�_�̓��̓C���s�[�_���X
//		double  VSWR_50[ ]  VSWR_75[ ]	�e���d�_��50����75���ɑ΂���VSWR
//--���ˊE�z��
//		complex TRAD[ ] FRAD[ ]			E�Ƃ�E�ӂ̕��ˊE
//		complex RRAD[ ] LRAD[ ]			ER ��EL �̕��ˊE
//		double TPHASE[ ] FPHASE[ ]		E�Ƃ�E�ӂ̈ʑ�
//		double RPHASE[ ] LPHASE[ ]		ER ��EL �̈ʑ�
//		double TGAI[ ] FGAI[ ]			E�Ƃ�E�ӂ̗���
//		double TFGAI[ ]					E�Ƃ�E�ӂ̍��v�̗���
//		double RLGAI[ ]					ER ��EL �̍��v�̗���
//		double RGAI[ ] LGAI[ ]			ER ��EL �̗���
//		double TPAT[ ] FPAT[ ]			E�Ƃ�E�ӂ̃p�^�[��
//		double RPAT[ ] LPAT[ ]			ER ��EL �̃p�^�[��
//		double AR[ ]					����
//--��{�����萔
//		double PI						��
//		double C						����
//		double e0						�^�󒆂̓�������0
//		double u0						�^�󒆂̗U�d����0
//		double R						1.0 �Œ�
//--���ȗ��p
//		complex J						�����P�ʁ@J=��-1
//--�ϕ��p
//		int GaussTenNo					�K�E�X�ϕ����_  [ 4�Œ�]�i�ʏ�p�j
//		int GaussTenSpe					�K�E�X�ϕ����_��[40�Œ�]�i���ٓ_�p�j
//--�v�Z�p
//		double k0						k0
//		complex uair					u(air)
//		double Beta;					��
//============================================================================
//============================================================================
