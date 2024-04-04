//	******************************************************************************************************
//	******************************************************************************************************
//	******																							******
//	******							 Numerical Mie calculation									    ******
//	******                           for Concentric Coated Sphere									******
//	******																							******
//	******************************************************************************************************
//	******************************************************************************************************
//	******																							******
//	******			A Comparison with "O.B.Toon & T.P. Ackerman, Algorithms for the calculation of	******
//	******			scattering by stratified spheres, Applied Optics 20(20), 1981, pp.3657-3660.	******
//	******																							******
//	******************************************************************************************************
//	******************************************************************************************************
//	******																							******
//	******						                           by Jianqi Shen   on Apr.10, 2006			******
//	******									               Email: jqshenk@163.com					******
//	******																							******
//	******************************************************************************************************
//	******************************************************************************************************

#include <cmath>
#include <cstdlib>
#include <cstdio>

const int j=100000;
double r_L11[j],i_L11[j],r_L12[j],i_L12[j],r_L22[j],i_L22[j],La1n[j];
double r_a[j],i_a[j],r_b[j],i_b[j],intens[3];

int MieAnBn_Coated(double alpha1,double rm1,double im1,double alpha2,double rm2,double im2,double k_coef[]);
void CX(double p,double q,double s,double t,double x[]);
void recursNFunc(int n,double rma,double ima,double r_Func,double i_Func,double y[]);
void recursXnz(int n,double rma,double ima,double r_Func,double i_Func,double y[]);
void FuncU(double r1,double i1,double rF1,double iF1,double r2,double i2,double rF2,double iF2,double y[]);
void FGL(double rma,double ima,int L,double ParaFGL[]);

void coated_mie(double alpha1, double rm1, double im1, double alpha2, double rm2, double im2, double k_coef[3])
{
	int n_max;
//	============================================================================================
	n_max=MieAnBn_Coated(alpha1,rm1,im1,alpha2,rm2,im2,k_coef);
//	===========================================================================================
}


//	********************************************************************************************************
//	********************************************************************************************************
//	*****************																		****************
//	*****************						MieCoefficients: an & bn						****************
//	*****************																		****************
//	*****************						   Kext,Ksca and Kabs							****************
//	*****************																		****************
//	********************************************************************************************************
//	********************************************************************************************************
int MieAnBn_Coated(double alpha1,double rm1,double im1,double alpha2,double rm2,double im2,double k_coef[])
{
	if (alpha1<alpha2)
	{
		double tmp;
		tmp=alpha1,alpha1=alpha2,alpha2=tmp;
		tmp=rm1,rm1=rm2,rm2=tmp;
		tmp=im1,im1=im2,im2=tmp;
	}

	int kim;
	if (im1<0||im2<0) im1=-fabs(im1),im2=-fabs(im2),kim=-1;
	else im1=fabs(im1),im2=fabs(im2),kim=1;
	int n=0,n_max=int(alpha1+7.5*pow(alpha1,0.34)+2.5);
	if (n_max>=j-1)
	{
		printf("\t\t\tCaution:\n\t\t\t\tThe parameter for the array is not large enough!!!");
		return (0);
	}
	double kext=0.0,ksca=0.0,com_x[4],com_y[2],ParaFGL[3];
	double alpha1s=pow(alpha1,2);
	double rm2a2=rm2*alpha2,rm1a2=rm1*alpha2,rm1a1=rm1*alpha1;
	double im2a2=im2*alpha2,im1a2=im1*alpha2,im1a1=im1*alpha1;
	double k0=0.0,k1=1.0,r1_tmp,i1_tmp,r2_tmp,i2_tmp,r3_tmp,i3_tmp;

//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//	[1] ==================================================================
	FGL(rm1a1,im1a1,n_max,ParaFGL);
	r_L11[n_max]=r2_tmp=ParaFGL[1],i_L11[n_max]=i2_tmp=ParaFGL[2];
	for (n=n_max;n>=1;n--)
	{
		CX(double(n),k0,rm1a1,im1a1,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(k1,k0,r1_tmp+r2_tmp,i1_tmp+i2_tmp,com_x);
		r_L11[n-1]=r2_tmp=r1_tmp-com_x[0],i_L11[n-1]=i2_tmp=i1_tmp-com_x[1];
	}

//	[2] ==================================================================
	FGL(rm1a2,im1a2,n_max,ParaFGL);
	r_L12[n_max]=r2_tmp=ParaFGL[1],i_L12[n_max]=i2_tmp=ParaFGL[2];
	for (n=n_max;n>=1;n--)
	{
		CX(double(n),k0,rm1a2,im1a2,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(k1,k0,r1_tmp+r2_tmp,i1_tmp+i2_tmp,com_x);
		r_L12[n-1]=r2_tmp=r1_tmp-com_x[0],i_L12[n-1]=i2_tmp=i1_tmp-com_x[1];
	}

//	[3] ==================================================================
	FGL(rm2a2,im2a2,n_max,ParaFGL);
	r_L22[n_max]=r2_tmp=ParaFGL[1],i_L22[n_max]=i2_tmp=ParaFGL[2];
	for (n=n_max;n>=1;n--)
	{
		CX(double(n),k0,rm2a2,im2a2,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(k1,k0,r1_tmp+r2_tmp,i1_tmp+i2_tmp,com_x);
		r_L22[n-1]=r2_tmp=r1_tmp-com_x[0],i_L22[n-1]=i2_tmp=i1_tmp-com_x[1];
	}

//	[4] ====================================================================
	FGL(alpha1,k0,n_max,ParaFGL);
	La1n[n_max]=ParaFGL[1];
	for (n=n_max;n>=1;n--) La1n[n-1]=double(n)/alpha1-1.0/(double(n)/alpha1+La1n[n]);

//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	double r_Ma1n,r_N12n,r_N11n,r_Na1n,r_Qn,i_Ma1n,i_N12n,i_N11n,i_Na1n,i_Qn;
//	M1(a1):
	r_Ma1n=pow(sin(alpha1)-alpha1*cos(alpha1),2)/(alpha1s+1.0);
	if (kim<0) i_Ma1n=((alpha1s-1.0)*sin(2.0*alpha1)*0.5+alpha1*cos(2.0*alpha1))/(alpha1s+1.0);
	else i_Ma1n=-((alpha1s-1.0)*sin(2.0*alpha1)*0.5+alpha1*cos(2.0*alpha1))/(alpha1s+1.0);

//	N0(a1),N0(m1a1),N0(m1a2):
	r_N12n=r_N11n=r_Na1n=0.0;
	if (kim<0) i_N12n=i_N11n=i_Na1n=-1.0;
	else i_N12n=i_N11n=i_Na1n=1.0;

//	Q0(m1a2m1a1):
	if (kim<=0)
	{
		CX(cos(2.0*rm1a2)-exp(2.0*im1a2),sin(2.0*rm1a2),cos(2.0*rm1a1)-exp(2.0*im1a1),sin(2.0*rm1a1),com_x);
		r_Qn=com_x[0]*exp(2.0*(im1a1-im1a2)),i_Qn=com_x[1]*exp(2.0*(im1a1-im1a2));
	}
	else
	{
		CX(exp(-2.0*im1a2)-cos(2.0*rm1a2),sin(2.0*rm1a2),exp(-2.0*im1a1)-cos(2.0*rm1a1),sin(2.0*rm1a1),com_x);
		r_Qn=com_x[0]*exp(2.0*(im1a2-im1a1)),i_Qn=com_x[1]*exp(2.0*(im1a2-im1a1));
	}

//	=======================================================================================================
//	=======================================================================================================
	for (n=1;n<=n_max;n++)
	{
		double nalpha=double(n)/alpha1;
//	[5] =========================================================================
		if (n>1)
		{
			CX(r_Ma1n,i_Ma1n,r_Na1n-nalpha,i_Na1n,com_x);
			r_Ma1n=com_x[0]*(La1n[n-1]-nalpha),i_Ma1n=com_x[1]*(La1n[n-1]-nalpha);
		}
//	[6] =========================================================================
		CX(k1,k0,nalpha-r_Na1n,-i_Na1n,com_x);
		r_Na1n=com_x[0]-nalpha,i_Na1n=com_x[1];
//	[7] =======================================================================
		recursNFunc(n,rm1a1,im1a1,r_N11n,i_N11n,com_y);
		r_N11n=com_y[0],i_N11n=com_y[1];
//	[8] =======================================================================
		recursNFunc(n,rm1a2,im1a2,r_N12n,i_N12n,com_y);
		r_N12n=com_y[0],i_N12n=com_y[1];
//	[9] =======================================================================
		recursXnz(n,rm1a1,im1a1,r_L11[n],i_L11[n],com_y);
		r1_tmp=com_y[0],i1_tmp=com_y[1];
		recursXnz(n,rm1a2,im1a2,r_N12n,i_N12n,com_y);
		r2_tmp=com_y[0],i2_tmp=com_y[1];
		CX(r1_tmp,i1_tmp,r2_tmp,i2_tmp,com_x);
		r3_tmp=com_x[2],i3_tmp=com_x[3];

		recursXnz(n,rm1a2,im1a2,r_L12[n],i_L12[n],com_y);
		r1_tmp=com_y[0],i1_tmp=com_y[1];
		recursXnz(n,rm1a1,im1a1,r_N11n,i_N11n,com_y);
		r2_tmp=com_y[0],i2_tmp=com_y[1];
		CX(r1_tmp,i1_tmp,r2_tmp,i2_tmp,com_x);
		r2_tmp=com_x[2],i2_tmp=com_x[3];

		CX(r3_tmp,i3_tmp,r2_tmp,i2_tmp,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(r_Qn,i_Qn,r1_tmp,i1_tmp,com_x);
		r_Qn=com_x[2],i_Qn=com_x[3];

//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//	[10] =====================================================================================
		double ru1a,ru2a,ru3a,ru4a,ru5a,ru6a,iu1a,iu2a,iu3a,iu4a,iu5a,iu6a;

		FuncU(rm2,im2,r_N12n,i_N12n,rm1,im1,r_L22[n],i_L22[n],com_y);
		ru1a=com_y[0],iu1a=com_y[1];

		FuncU(1.0,0.0,r_L11[n],i_L11[n],rm1,im1,La1n[n],0.0,com_y);
		ru2a=com_y[0],iu2a=com_y[1];

		FuncU(rm2,im2,r_L12[n],i_L12[n],rm1,im1,r_L22[n],i_L22[n],com_y);
		ru3a=com_y[0],iu3a=com_y[1];

		FuncU(1.0,0.0,r_N11n,i_N11n,rm1,im1,La1n[n],0.0,com_y);
		ru4a=com_y[0],iu4a=com_y[1];

		FuncU(1.0,0.0,r_L11[n],i_L11[n],rm1,im1,r_Na1n,i_Na1n,com_y);
		ru5a=com_y[0],iu5a=com_y[1];

		FuncU(1.0,0.0,r_N11n,i_N11n,rm1,im1,r_Na1n,i_Na1n,com_y);
		ru6a=com_y[0],iu6a=com_y[1];

//	[11] =====================================================================================
		double ru1b,ru2b,ru3b,ru4b,ru5b,ru6b,iu1b,iu2b,iu3b,iu4b,iu5b,iu6b;
		FuncU(rm1,im1,r_N12n,i_N12n,rm2,im2,r_L22[n],i_L22[n],com_y);
		ru1b=com_y[0],iu1b=com_y[1];

		FuncU(rm1,im1,r_L11[n],i_L11[n],1.0,0.0,La1n[n],0.0,com_y);
		ru2b=com_y[0],iu2b=com_y[1];

		FuncU(rm1,im1,r_L12[n],i_L12[n],rm2,im2,r_L22[n],i_L22[n],com_y);
		ru3b=com_y[0],iu3b=com_y[1];

		FuncU(rm1,im1,r_N11n,i_N11n,1.0,0.0,La1n[n],0.0,com_y);
		ru4b=com_y[0],iu4b=com_y[1];

		FuncU(rm1,im1,r_L11[n],i_L11[n],1.0,0.0,r_Na1n,i_Na1n,com_y);
		ru5b=com_y[0],iu5b=com_y[1];

		FuncU(rm1,im1,r_N11n,i_N11n,1.0,0.0,r_Na1n,i_Na1n,com_y);
		ru6b=com_y[0],iu6b=com_y[1];

//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//	[12] =========================================================================================
		CX(r_Qn,i_Qn,ru3a,iu3a,com_x);
		r1_tmp=com_x[2],i1_tmp=com_x[3];
		FuncU(ru1a,iu1a,ru2a,iu2a,r1_tmp,i1_tmp,ru4a,iu4a,com_y);
		r2_tmp=com_y[0],i2_tmp=com_y[1];
		FuncU(ru1a,iu1a,ru5a,iu5a,r1_tmp,i1_tmp,ru6a,iu6a,com_y);
		r3_tmp=com_y[0],i3_tmp=com_y[1];
		CX(r2_tmp,i2_tmp,r3_tmp,i3_tmp,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(r_Ma1n,i_Ma1n,r1_tmp,i1_tmp,com_x);
		r_a[n]=com_x[2],i_a[n]=com_x[3];

//	[13] =========================================================================================
		CX(r_Qn,i_Qn,ru3b,iu3b,com_x);
		r1_tmp=com_x[2],i1_tmp=com_x[3];
		FuncU(ru1b,iu1b,ru2b,iu2b,r1_tmp,i1_tmp,ru4b,iu4b,com_y);
		r2_tmp=com_y[0],i2_tmp=com_y[1];
		FuncU(ru1b,iu1b,ru5b,iu5b,r1_tmp,i1_tmp,ru6b,iu6b,com_y);
		r3_tmp=com_y[0],i3_tmp=com_y[1];
		CX(r2_tmp,i2_tmp,r3_tmp,i3_tmp,com_x);
		r1_tmp=com_x[0],i1_tmp=com_x[1];
		CX(r_Ma1n,i_Ma1n,r1_tmp,i1_tmp,com_x);
		r_b[n]=com_x[2],i_b[n]=com_x[3];

//	[14] =============================================================================
		double para=2.0*double(n)+1.0;
		kext+=para*(r_a[n]+r_b[n]);
		ksca+=para*(pow(r_a[n],2)+pow(i_a[n],2)+pow(r_b[n],2)+pow(i_b[n],2));

	}
	if ((rm1==1.0)&&(im1==0.0)) alpha1s=pow(alpha2,2);
	k_coef[0]=kext*2.0/alpha1s,k_coef[1]=ksca*2.0/alpha1s,k_coef[2]=k_coef[0]-k_coef[1];
	if (fabs(k_coef[2]/k_coef[0])<1e-10) k_coef[2]=0.0;
	return (n_max);
}

//	********************************************************************************************************
//	********************************************************************************************************
//	*****************																		****************
//	*****************							Complex Calculation							****************
//	*****************																		****************
//	********************************************************************************************************
//	********************************************************************************************************
void CX(double p,double q,double s,double t,double x[])
{
	double modl=s*s+t*t;
	x[0]=(p*s+q*t)/modl,x[1]=(q*s-p*t)/modl,x[2]=p*s-q*t,x[3]=p*t+q*s;
}


//	********************************************************************************************************
//	********************************************************************************************************
//	********************************************************************************************************
//	********************************************************************************************************
void recursNFunc(int n,double rma,double ima,double r_Func,double i_Func,double y[])
{
	double rtmp,itmp,k1=1.0,k0=0.0,com_x[4];
	CX(double(n),k0,rma,ima,com_x);
	rtmp=com_x[0],itmp=com_x[1];
	CX(k1,k0,rtmp-r_Func,itmp-i_Func,com_x);
	y[0]=com_x[0]-rtmp,y[1]=com_x[1]-itmp;
	return;
}

//	********************************************************************************************************
//	********************************************************************************************************
//	*****************																		****************
//	*****************																		****************
//	********************************************************************************************************
//	********************************************************************************************************
void recursXnz(int n,double rma,double ima,double r_Func,double i_Func,double y[])
{
	double com_x[4];
	CX(n,0.0,rma,ima,com_x);
	y[0]=com_x[0]+r_Func,y[1]=com_x[1]+i_Func;
	return;
}

//	********************************************************************************************************
//	********************************************************************************************************
//	*****************																		****************
//	*****************																		****************
//	********************************************************************************************************
//	********************************************************************************************************
void FuncU(double r1,double i1,double rF1,double iF1,double r2,double i2,double rF2,double iF2,double y[])
{
	double com_x[4];
	CX(r1,i1,rF1,iF1,com_x);
	y[0]=com_x[2],y[1]=com_x[3];
	CX(r2,i2,rF2,iF2,com_x);
	y[0]-=com_x[2],y[1]-=com_x[3];
	return;
}


//	********************************************************************************************************
//	********************************************************************************************************
//	*****************																		****************
//	*****************						Lentz: L(n_max)									****************
//	*****************																		****************
//	********************************************************************************************************
//	********************************************************************************************************
void FGL(double rma,double ima,int L,double ParaFGL[])
{
	double rfl3,ifl3,com_x[4];
	double k0=0.0,k1=2.0*L+1.0,k2=-2.0*L-3.0;
	CX(k1,k0,rma,ima,com_x);
	double rwl1=com_x[0],iwl1=com_x[1];
	CX(k2,k0,rma,ima,com_x);
	double rwl2=com_x[0],iwl2=com_x[1];
	double rfl2=rwl1+rma/k2,ifl2=iwl1+ima/k2;
	double ro2=rwl2+rma/k1,io2=iwl2+ima/k1;
	double rp2=rwl2,ip2=iwl2;
	double rwlk=rwl2,iwlk=iwl2;
	for (int k=3;;k++)
	{
		double parak=-(2.0*(L+k)-1)/(2.0*(L+k)-3);
        rwlk=parak*rwlk,iwlk=parak*iwlk;
        k1=1.0;
        CX(k1,k0,ro2,io2,com_x);
        double ro3=rwlk+com_x[0],io3=iwlk+com_x[1];
        CX(k1,k0,rp2,ip2,com_x);
        double rp3=rwlk+com_x[0],ip3=iwlk+com_x[1];
        CX(ro3,io3,rp3,ip3,com_x);
		double rp4=com_x[0],ip4=com_x[1];
        CX(rfl2,ifl2,rp4,ip4,com_x);
		rfl3=com_x[2],ifl3=com_x[3];
        if (fabs(1-sqrt(rp4*rp4+ip4*ip4))<1e-20) break;
		else ro2=ro3,io2=io3,rp2=rp3,ip2=ip3,rfl2=rfl3,ifl2=ifl3;
	}
	k1=-L,k0=0.0;
	CX(k1,k0,rma,ima,com_x);
	ParaFGL[1]=rfl3+com_x[0],ParaFGL[2]=ifl3+com_x[1];
}