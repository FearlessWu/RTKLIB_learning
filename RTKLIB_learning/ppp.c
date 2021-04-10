/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2016 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS ���ʵ�����ת���� tide model
*           -DOUTSTAT_AMB output ambiguity ģ�� parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Co
entions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solutin problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[]="$Id:$";

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations ���������� */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution ������� */
#define MIN_NSAT_SOL 4              /* min satellite number for solution �������������������Ŀ */
#define THRES_REJECT 4.0            /* reject threshold����ֵ��of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance���仯�� receiver position (m^2)��ʼ�仯���ջ���λ�� */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) ��ʼ�仯���ջ�ʱ�� */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2)��ʼ�仯ZTD��ֵ */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2)�ݶ� gradient���ݶȣ�*/
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2)��λƫ�� */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay ������ӳ�*/
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity��ʪ�ȣ� for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ���� ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular ��ǵ� velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

/* number and index of states ������״ָ̬��  �궨��ĺ���*/
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)	//Ƶ��������ȡ
#define NP(opt)     ((opt)->dynamics?9:3)						//��̬ģ��
#define NC(opt)     (NSYS)										//����ϵͳ������ȡ
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))//������ģ��
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)		//������Ƿ���� ����ģ��
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))                           //9:3+NSYS
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)//3+1+(0/1/3)+i+1
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))//3+1+(0/1/3)+maxsat
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

/* standard deviation of state /����� �̶���ı�׼�� -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);//����ǹ̶�״̬��ô������  nx--������   SOLQ_FIX��ʾ1Ϊ�̶�״̬  sol.stat--�������  Pa�̶���Э����
    return SQRT(rtk->P[i+i*rtk->nx]);// �����Э����--P
}
/* write solution status for PPP /���ܶ�λ�����д��---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)//buff����
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],*x;//
    int i,j,week;
    char id[32],*p=buff;
    
    if (!rtk->sol.stat) return 0;//sol_stat  �������Ϊ��ʱ 
    
    trace(3,"pppoutstat:\n");// void trace(int leve,const char *format,....)?????
    
    tow=time2gpst(rtk->sol.time,&week);// void time2gpst(gtime_t ,char ,int)
    
    x=rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;//���̶���x=�̶��⣬��֮x=�����
    
    /* receiver position */
    p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2));
    
    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {            //opt.dynamics(0:none,1:velociy,2:accel)������0��--����ѡ���̬��ʽ
        ecef2pos(rtk->sol.rr,pos);      //sol.rr��ά���������������
        ecef2enu(pos,rtk->x+3,vel);    //rtk->x �����
        ecef2enu(pos,rtk->x+6,acc);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                   "%.4f,%.5f,%.5f,%.5f\n",week,tow,rtk->sol.stat,vel[0],vel[1],
                   vel[2],acc[0],acc[1],acc[2],0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks  ���ջ�ʱ��*/

    i=IC(0,&rtk->opt);
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,x[i]*1E9/CLIGHT,x[i+1]*1E9/CLIGHT,
               STD(rtk,i)*1E9/CLIGHT,STD(rtk,i+1)*1E9/CLIGHT);
    
    /* tropospheric parameters  ��������� */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {  //opt.tropopt--������ѡ�TROPOPT_EST=3 ����troposphere option: ZTD estimation
        i=IT(&rtk->opt);                                                  // TROPOPT_ESTG 4    troposphere option: ZTD+grad estimation 
        p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                   1,x[i],STD(rtk,i));
    }
    if (rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n",week,tow,
                   rtk->sol.stat,1,x[i+1],x[i+2],STD(rtk,i+1),STD(rtk,i+2));  /* double std [MAXSAT][3]-- fcb std-dev (cyc) */
    }
    /* ionosphere parameters �������� */
    if (rtk->opt.ionoopt==IONOOPT_EST) {   //#define IONOOPT_EST 4   ionosphere option: estimation 
        for (i=0;i<MAXSAT;i++) {        //MAXSAT--���������Ŀ
            ssat=rtk->ssat+i;            /* satellite status /����״̬  ������λ�Ǻ͸߶�*/
            if (!ssat->vs) continue;       /*ssat->vs=0 ���� valid satellite flag single /��Ч���Ǳ�־�� */
            j=II(i+1,&rtk->opt);
            if (rtk->x[j]==0.0) continue;
            satno2id(i+1,id);
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,rtk->ssat[i].azel[0]*R2D,
                       rtk->ssat[i].azel[1]*R2D,x[j],STD(rtk,j));         /*FLOAT azel[2]--��λ�Ǻ͸߶� azimuth/elevation (rad) */
        }
    }
#ifdef OUTSTAT_AMB       /* ambiguity parameters ģ���Ȳ���*/
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NF(&rtk->opt);j++) {
        k=IB(i+1,j,&rtk->opt);
        if (rtk->x[k]==0.0) continue;
        satno2id(i+1,id);// void satno2id(int sat, char *id);
        p+=sprintf(p,"$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,id,j+1,x[k],STD(rtk,k));//?????
    }
#endif
    return (int)(p-buff);//??????
}
/* exclude meas of eclipsing satellite (block IIA) /�ų���ʧ/������/ʳ���� ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");// void trace(int leve,const char *format,....)?????
    
    /* unit vector of sun direction (ecef)̫������ĵ�λʸ�� */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);//���̫�����������ڵع�����ϵ��λ������
	//obs[0].time--���ջ�����ʱ�䣻erpv--������ת������rsun--̫���ع���������
	/*void sunmoonpos(gtime_t tutc,const double *erpv, double *rsun,double *rmoon,double *gmst);   */

    normv3(rsun,esun);//����̫���ķ�������
	// int  normv3(const double *a, double *b)
    
    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;//antenna type
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;//���벻С����
        
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;//���ַ���str1�в����Ƿ����ַ���str2�� ����У���str1�е�str2λ���𣬷���str1��ָ�룬���û�У�����null
        
        /* sun-earth-satellite angle ̫��-����-���ǽǶ�*/
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);//cosa��Χ��-1��1֮�䣬�ϻ�
        ang=acos(cosa);//�������-̫���Ƕȣ��������ǶȲ�����pi/2��
        
        /* test eclipse  �ص� */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;//���ǵ�������Ϊ��
    }
}
/* nominal yaw-angle /��Сƫ�Ǵ��� ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite /����ƫ�Ǵ��� return yaw ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model /������̬ģ��  --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;
    
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);
    
    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;
    
    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;
    
    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model ��λ����ģ�ͣ�����-------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;
    
    if (opt<=0) return 1; /* no phase windup */
    
    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;
    
    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;
    
    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */
    
    /* phase windup effect */
    cross3(ek,eys,eks);			
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;
    
    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* measurement error variance ��������------------------------------------------------*/
static double varerr(int sat, int sys, double el, int freq, int type,
                     const prcopt_t *opt)
{
    double fact=1.0,sinel=sin(el);
    
    if (type==1) fact*=opt->eratio[freq==0?0:1];
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    
    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (freq==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
    return SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
}
/* initialize state and covariance ��ʼ״̬��Э����-------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

/*x[]--����⣻nx--����������δ֪��������*/
/* geometry-free phase measurement �޼�����λ����-------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0) return 0.0;
    return lam[0]*obs->L[0]-lam[i]*obs->L[i];
}
/* Melbourne-Wubbena linear combination ������� --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0||
        obs->P[0]==0.0||obs->P[i]==0.0) return 0.0;
    return lam[0]*lam[i]*(obs->L[0]-obs->L[i])/(lam[i]-lam[0])-
           (lam[i]*obs->P[0]+lam[0]*obs->P[i])/(lam[i]+lam[0]);
}
/* antenna corrected measurements  ���߸������� --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    const double *lam=nav->lam[obs->sat-1];
    double C1,C2;
    int i,sys;
    
    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        if (lam[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) continue;
        
        /* antenna phase center and phase windup correction ������λ���ĺ���λ�������*/
        L[i]=obs->L[i]*lam[i]-dants[i]-dantr[i]-phw*lam[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];
        
        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) �����ƫ����� */
        if (obs->code[i]==CODE_L1C) {
            P[i]+=nav->cbias[obs->sat-1][1];
        }
        else if (obs->code[i]==CODE_L2C||obs->code[i]==CODE_L2X||
                 obs->code[i]==CODE_L2L||obs->code[i]==CODE_L2S) {
            P[i]+=nav->cbias[obs->sat-1][2];
#if 0
            L[i]-=0.25*lam[i]; /* 1/4 cycle-shift */
#endif
        }
    }
    /* iono-free LC �޵����ģ��*/
    *Lc=*Pc=0.0;
    sys=satsys(obs->sat,NULL);
    i=(sys&(SYS_GAL|SYS_SBS))?2:1; /* L1/L2 or L1/L5 */
    if (lam[0]==0.0||lam[i]==0.0) return;
    
    C1= SQR(lam[i])/(SQR(lam[i])-SQR(lam[0]));
    C2=-SQR(lam[0])/(SQR(lam[i])-SQR(lam[0]));
    
#if 0
    /* P1-P2 dcb correction (P1->Pc,P2->Pc) */
    if (sys&(SYS_GPS|SYS_GLO|SYS_QZS)) {
        if (P[0]!=0.0) P[0]-=C2*nav->cbias[obs->sat-1][0];
        if (P[1]!=0.0) P[1]+=C1*nav->cbias[obs->sat-1][0];
    }
#endif
    if (L[0]!=0.0&&L[i]!=0.0) *Lc=C1*L[0]+C2*L[i];
    if (P[0]!=0.0&&P[i]!=0.0) *Pc=C1*P[0]+C2*P[i];
}
/* detect cycle slip by LLI ������� --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j;
    
    trace(3,"detslp_ll: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        
        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);
        
        rtk->ssat[obs[i].sat-1].slip[j]=1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j;
    
    trace(3,"detslp_gf: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        if ((g1=gfmeas(obs+i,nav))==0.0) continue;
        
        g0=rtk->ssat[obs[i].sat-1].gf;
        rtk->ssat[obs[i].sat-1].gf=g1;
        
        trace(4,"detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n",obs[i].sat,g0,g1);
        
        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            trace(3,"detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",
                  obs[i].sat,g0,g1);
            
            for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double w0,w1;
    int i,j;
    
    trace(3,"detslp_mw: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((w1=mwmeas(obs+i,nav))==0.0) continue;
        
        w0=rtk->ssat[obs[i].sat-1].mw;
        rtk->ssat[obs[i].sat-1].mw=w1;
        
        trace(4,"detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f\n",obs[i].sat,w0,w1);
        
        if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
            trace(3,"detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n",
                  obs[i].sat,w0,w1);
            
            for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* temporal update of position �ݴ�λ�õĸ���-----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    int i;
    
    trace(3,"udpos_ppp:\n");
    
    /* fixed mode �̶�ģ��*/
    if (rtk->opt.mode==PMODE_PPP_FIXED) {             // when the mode is ppp_fixedn
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);

		//������1E-8ΪЭ����Խ��󣬾����СΪ3*nx

        return;
    }//this mode is not been done yet.
	
	/* initialize position for first epoch ��һ��Ԫλ�õĳ�ʼ��*/
    if (norm(rtk->x,3)<=0.0) {    //norm function :����Ԫ��ƽ���Ϳ������˴�Ϊ����С��0��
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);//receiver position vector is 60^2
    }


    /* static ppp mode��̬pppģ�� */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {//ģ��Ϊ��̬ppp��ʱ��
        for (i=0;i<3;i++) {
            rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);//ԭ����Э����+����*��Ԫʱ���
        }
        return;
    }
	//��ģ��Ϊ��̬pppʱ�����ɵ�Э������ΪP[i,i]=��������^2*|delta t|


    /* kinmatic mode dynamics */
    for (i=0;i<3;i++) {
        initx(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
}
//receiver position vector is 60^2






/* temporal update of clock ʱ�ӵ�ʱ�����--------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
    double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise)������ ÿһ��Ԫ�ӳ�ʼ��  */
    for (i=0;i<NSYS;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {  //����������Ϊ��������
            /* time of prec ephemeris is based gpst ���������ǻ���GPSʱ��*/
            /* negelect receiver inter-system bias ���Խ��ջ����ϵͳ��� */
            dtr=rtk->sol.dtr[3];//���ջ��Ӳ�
			//dtr = rtk->sol.dtr[0];//���ջ��Ӳ�
        }
        else {
            dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];//�Ǿ�������
        }
        initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));//�Ӳ����ʼ����Э����
//�Ӳ���Ϊc*dtr,Э����Ϊ60��ƽ����
    }
}
/* temporal update of tropospheric parameters �����������ʱ�����--------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&rtk->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    
    if (rtk->x[i]==0.0) {//�������㸡���Ϊ��
        ecef2pos(rtk->sol.rr,pos);//�ع�����ת��������꣬��������ʽΪ��latitude��longtitude��height��
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);//EGNOSģ�Ͷ������ӳٸ���
        initx(rtk,ztd,var,i);//�����㸡����Э�����ʼ���������Ϊ�������ӳٸ�������Э����Ϊ����������,�����СΪ1*1
        
		if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(rtk,1E-6,VAR_GRA,j);
			//����ģʽ֮�£��ټ�����両���Ϊ1E-6��Э�����Ϊ0.01^2,�����СΪ2*2
        }
    }
    else {
        rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);//prn[2]--������������׼ƫ��ֵ���õ�Э����Ϊ��org+��std*t��^2
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);//prn[2]--������������׼ƫ��ֵ���õ�Э����Ϊ��org+��std*0.1*t��^2
            }
        }
    }
}
/* temporal update of ionospheric parameters �����Ĳ���ʱ����� ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double ion,sinel,pos[3],*azel;//lam�ز�������pos���ջ�λ�á�azel�����춥��
    char *p;
    int i,j,k,gap_resion=GAP_RESION;//���õ�������120
    
    trace(3,"udiono_ppp:\n");
    
    if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) {//strstr()�������ҳ������ַ���ͬ��λ�ã��ҵ����ظ�λ��ָ�룬���򷵻��㡣
        sscanf(p,"-GAP_RESION=%d",&gap_resion);
    }
    for (i=0;i<MAXSAT;i++) {
        j=II(i+1,&rtk->opt);
        if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) {//����ⲻ�����㣬����λ�жϴ�������120��������ʼΪ0
            rtk->x[j]=0.0;
        }
    }
    for (i=0;i<n;i++) {
        j=II(obs[i].sat,&rtk->opt);
        if (rtk->x[j]==0.0) {
            k=satsys(obs[i].sat,NULL)==SYS_GAL?2:1;
            lam=nav->lam[obs[i].sat-1];//�ز�����
            if (obs[i].P[0]==0.0||obs[i].P[k]==0.0||lam[0]==0.0||lam[k]==0.0) {//PΪ α��۲�����
                continue;
            }
            ion=(obs[i].P[0]-obs[i].P[k])/(1.0-SQR(lam[k]/lam[0]));
            ecef2pos(rtk->sol.rr,pos);//�ع�����ת�������
            azel=rtk->ssat[obs[i].sat-1].azel;//��λ������
            ion/=ionmapf(pos,azel);//����㴩͸�㣨IPP��λ�ú���б����
            initx(rtk,ion,VAR_IONO,j);
        }
        else {
            sinel=sin(MAX(rtk->ssat[obs[i].sat-1].azel[1],5.0*D2R));
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]/sinel)*fabs(rtk->tt);
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
    int i=ID(&rtk->opt);
    
    trace(3,"uddcb_ppp:\n");
    
    if (rtk->x[i]==0.0) {
        initx(rtk,1E-6,VAR_DCB,i);//�����Ϊ1E-6��Э����Ϊ30ƽ����
    }
}
/* temporal update of phase biases ��λƫ��-------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double ion,dantr[NFREQ]={0},dants[NFREQ]={0};
    int i,j,k,l,f,sat,slip[MAXOBS]={0},clk_jump=0;
    
    trace(3,"udbias  : n=%d\n",n);
    
    /* handle day-boundary clock jump//�����ձ߽�ʱ����ת */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;//��ʼ��������ʶ
    }
    /* detect cycle slip by LLI ����̽�� */
    detslp_ll(rtk,obs,n);
    
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk,obs,n,nav);
    
    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(rtk,obs,n,nav);
//������� ģ�ͣ��������������rtk->ssat[obs[i].sat-1].slip[j]|=1
    
    ecef2pos(rtk->sol.rr,pos);
    
    for (f=0;f<NF(&rtk->opt);f++) {
        
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));//����۲����ݵ��жϼ�����ʧЧ�����ø������Э������
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);
            corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,dantr,dants,
                      0.0,L,P,&Lc,&Pc);//���߸���������������L��P�Ĺ۲�ֵ��
            
            bias[i]=0.0;
            
            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc-Pc;
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];//�޵������ϵ�ѡ��
            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                slip[i]=rtk->ssat[sat-1].slip[f];
                l=satsys(sat,NULL)==SYS_GAL?2:1;
                lam=nav->lam[sat-1];
                if (obs[i].P[0]==0.0||obs[i].P[l]==0.0||
                    lam[0]==0.0||lam[l]==0.0||lam[f]==0.0) continue;
                ion=(obs[i].P[0]-obs[i].P[l])/(1.0-SQR(lam[l]/lam[0]));
                bias[i]=L[f]-P[f]+2.0*ion*SQR(lam[f]/lam[0]);
            }
            if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue;
            
            offset+=bias[i]-rtk->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=IB(i+1,f,&rtk->opt);
                if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);
            
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
            
            if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;
            
            /* reinitialize phase-bias if detecting cycle slip */
            initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt));
            
            /* reset fix flags */
            for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;
            
            trace(5,"udbias_ppp: sat=%2d bias=%.3f\n",sat,bias[i]);
        }
    }
}
/* temporal update of states ״̬�ĸ���--------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(3,"udstate_ppp: n=%d\n",n);
    
    /* temporal update of positionλ�ø��� */
    udpos_ppp(rtk);
	
/*//λ�ø��°������ջ���ά��������븡����Э������ĸ���
when the mode is ppp_fixed��������1E-8ΪЭ����Խ��󣬾����СΪ3*nx,��һ��Ԫλ�õĳ�ʼ��,
 receiver position vector is 60^2;
ģ��Ϊ��̬ppp��ʱ,���ɵ�Э������ΪP[i,i]=��������^2*|delta t|;
��̬ʱ��receiver position vector is 60^2 */
    

    /* temporal update of clock �Ӳ����*/
    udclk_ppp(rtk);
    /*�Ӳ���£��Ӳ���Ϊc*dtr,Э����Ϊ60��ƽ���������СΪ[NSYS*NSYS]*/

    /* temporal update of tropospheric parameters ������������� */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
	/*��������£�
	��������ѡ��ΪEstimate ZTD����Estimate ZTD+Grad
	�����Ϊ�������ӳٸ�����ztd��EGNOSģ�ͣ���
	ԭ�������Ϊ0ʱ����x[i]=0:
	                         Э����Ϊ����������var,�����СΪ1*1
	 ������ѡ��ΪTROPOPT_ESTGʱ���ټ�����両���Ϊ1E-6��Э�����Ϊ0.01^2,�����СΪ2*2
	else��
	     �õ�Э����Ϊ��org+��std*t��^2
		 ������ѡ��ΪTROPOPT_ESTGʱ��Э�������ټ�����õ�Э����Ϊ��org+��std*0.0*t��^2
	*/



    /* temporal update of ionospheric parameters ������������ */
    if (rtk->opt.ionoopt==IONOOPT_EST) {
        udiono_ppp(rtk,obs,n,nav);
    }


    /* temporal update of L5-receiver-dcb parameters���ƫ����� */
    if (rtk->opt.nf>=3) {  //��Ƶ��������3��ʱ
        uddcb_ppp(rtk);//�����Ϊ1E-6��Э����Ϊ30ƽ��
    }


    /* temporal update of phase-bias //��λƫ�����*/
    udbias_ppp(rtk,obs,n,nav);//ģ���ȵĸ����=bais[i]�Լ�Э����60square

}

/* satellite antenna phase center variation ����������λ���ı���----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;
    
    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;//euΪ���Ƿ������ң�ezΪ���ջ��Ե���Ϊԭ��ķ�������
    
    cosa=dot(eu,ez,3);//
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);//�Ƕ�
    
    antmodel_s(pcv,nadir,dant);//���߲�������dant
}
/* precise tropospheric model ��ȷ������ģ��------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;
    
    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);
    
    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);
    
    if (azel[1]>0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    double trp[3]={0},std[3];
    
    if (opt->tropopt==TROPOPT_SAAS) {//��˹��Ī��ģ��
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);//dtrpΪ��������ӳ�
        *var=SQR(ERR_SAAS);//varΪ����
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {//SBASģ��
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {//ZTD/ZTD+GRAD����
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);//���ܶ�����ģ��
        return 1;
    }
    if (opt->tropopt==TROPOPT_ZTD) {
        if (pppcorr_trop(&nav->pppcorr,time,pos,trp,std)) {
            *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);//�����������
            *var=SQR(dtdx[0]*std[0]);
            return 1;
        }
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    static double iono_p[MAXSAT]={0},std_p[MAXSAT]={0};
    static gtime_t time_p;
    
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_EST) {
        *dion=x[II(sat,opt)];
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_STEC) {
        if (timediff(time,time_p)!=0.0&&
            !pppcorr_stec(&nav->pppcorr,time,pos,iono_p,std_p)) return 0;
        if (iono_p[sat-1]==0.0||std_p[sat-1]>0.1) return 0;
        time_p=time;
        *dion=iono_p[sat-1];
        *var=SQR(std_p[sat-1]);
        return 1;
    }
    return 0;
}
/* constraint to local correction �ֲ�����Լ��--------------------------------------------*/
static int const_corr(const obsd_t *obs, int n, const int *exc,
                      const nav_t *nav, const double *x, const double *pos,
                      const double *azel, rtk_t *rtk, double *v, double *H,
                      double *var)
{
    gtime_t time=obs[0].time;
    double trop[3],std_trop[3],iono[MAXSAT],std_iono[MAXSAT];
    int i,j,k,sat,nv=0;
    
    /* constraint to external troposphere correction */
    if ((rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)&&
        pppcorr_trop(&nav->pppcorr,time,pos,trop,std_trop)) {//���������
        
        for (i=0;i<(rtk->opt.tropopt==TROPOPT_EST?1:3);i++) {
            if (std_trop[i]==0.0) continue;
            j=IT(&rtk->opt)+i;
            v[nv]=trop[i]-x[j];//�����������
            for (k=0;k<rtk->nx;k++) H[k+nv*rtk->nx]=k==j?1.0:0.0;
            var[nv++]=SQR(std_trop[i]);//�����㷽��
        }
    }
    /* constraint to external ionosphere correction */
    if (rtk->opt.ionoopt==IONOOPT_EST&&
        pppcorr_stec(&nav->pppcorr,time,pos,iono,std_iono)) {
        
        for (i=0;i<n;i++) {
            sat=obs[i].sat;
            if (exc[i]||iono[sat-1]==0.0||std_iono[sat-1]>0.5) continue;
            j=II(sat,&rtk->opt);
            v[nv]=iono[sat-1]-x[j];
            for (k=0;k<rtk->nx;k++) H[k+nv*rtk->nx]=k==j?1.0:0.0;
            var[nv++]=SQR(std_iono[sat-1]);
        }
    }
    return nv;
}
/* phase and code residuals ��λ�ͱ���ʣ�����--------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel)//rs--����λ�����ٶȣ�dts--�����Ӳvar_rs--����λ���Ӳ����svh--���ǽ���״����dr--��ϫ���������
	//exc--�����Ƿ��ų���x--δ֪������ֵ��v--������H--��ƾ����ת�ã�R--��������
{
    const double *lam;//����
    prcopt_t *opt=&rtk->opt;
    double y,r,cdtr,bias,C,rr[3],pos[3],e[3],dtdx[3],L[NFREQ],P[NFREQ],Lc,Pc;
    double var[MAXOBS*2],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,dcb;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    double ve[MAXOBS*2*NFREQ]={0},vmax=0;
    char str[32];
    int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
    int i,j,k,sat,sys,nv=0,nx=rtk->nx,stat=1;
    
    time2str(obs[0].time,str,2);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;//��ÿ�����ǽ�����Ч���Ǳ�־���г�ʼ������ʼ��ֵΪ0��Ƶ�������3��
    
    for (i=0;i<3;i++) rr[i]=x[i]+dr[i];//drΪ��ϫ��������
    ecef2pos(rr,pos);//���ջ��ع�����rrת���������ϵpos
    
    for (i=0;i<n&&i<MAXOBS;i++) {//��ÿ�����ǽ��в���
        sat=obs[i].sat;//satellite number
        lam=nav->lam[sat-1];//�ز�����
        if (lam[j/2]==0.0||lam[0]==0.0) continue;//j=2��������Ϊ����������ѭ��
        
        if ((r=geodist(rs+i*6,rr,e))<=0.0||//rs+i*6��ʾ�����i*6+1��Ԫ�أ��˺�����ý��ջ������ǵķ������ң�rΪ���ϵ�����ת�����ĺ�ľ���
            satazel(pos,e,azel+i*2)<opt->elmin) {//������ǵķ�λ�Ǻ͸߶Ƚ�azel��e���ع������µĽ��ջ�����
            exc[i]=1;//�ų�����
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||//1.ת�����ǺŻ�ȡ����ϵͳ 2.��Ч���Ǳ�־
            satexclude(obs[i].sat,svh[i],opt)||exc[i]) {//3.�ų��߶Ƚǲ����ϵ�����
            exc[i]=1;
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&vart)||//��ö������������Э����
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
            continue;
        }
        /* satellite and receiver antenna model */
        if (opt->posopt[0]) satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);//dants�������߲�������
        antmodel(opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr);//dantr���ջ����߸�����
        
        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {//��λ���Ƹ���
            continue;
        }
        /* corrected phase and code measurements ��λ�����������*/
        corr_meas(obs+i,nav,azel+i*2,&rtk->opt,dantr,dants,
                  rtk->ssat[sat-1].phw,L,P,&Lc,&Pc);
        
        /* stack phase and code residuals {L1,P1,L2,P2,...}������λ�����в� */
        for (j=0;j<2*NF(opt);j++) {
            
            dcb=bias=0.0;
            
            if (opt->ionoopt==IONOOPT_IFLC) {//ionosphere option: L1/L2 or L1/L5 iono-free LC 
                if ((y=j%2==0?Lc:Pc)==0.0) continue;//j=0,2Ϊ�ز���λ
            }
            else {
                if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;
                
                /* receiver DCB correction for P2 */
                if (j/2==1) dcb=-nav->rbias[0][sys==SYS_GLO?1:0][0];//rbias[64][2][3];���ջ���ƫ��
            }
            C=SQR(lam[j/2]/lam[0])*ionmapf(pos,azel+i*2)*(j%2==0?-1.0:1.0);//pos���ջ��������
            for (k=0;k<nx;k++) H[k+nx*nv]=k<3?-e[k]:0.0;//��ǰ������ϵ������Ϊ���������ĵ��෴��
            
            /* receiver clock */
            k=sys==SYS_GLO?1:(sys==SYS_GAL?2:(sys==SYS_CMP?3:0));
            cdtr=x[IC(k,opt)];//���ջ��Ӳ��ֵ
            H[IC(k,opt)+nx*nv]=1.0;//�Ӳ�ϵ������
            
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    H[IT(opt)+k+nx*nv]=dtdx[k];
                }
            }
            if (opt->ionoopt==IONOOPT_EST) {
                if (rtk->x[II(sat,opt)]==0.0) continue;
                H[II(sat,opt)+nx*nv]=C;
            }
            if (j/2==2&&j%2==1) { /* L5-receiver-dcb */
                dcb+=rtk->x[ID(opt)];
                H[ID(opt)+nx*nv]=1.0;
            }
            if (j%2==0) { /* phase bias */
                if ((bias=x[IB(sat,j/2,opt)])==0.0) continue;
                H[IB(sat,j/2,opt)+nx*nv]=1.0;//IB��λ���ڱ䣬��֤��ģ�����ھ����λ��
            }
            /* residual */
            v[nv]=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias);//��λ��α��Ĳв�
            
            if (j%2==0) rtk->ssat[sat-1].resc[j/2]=v[nv];//��λ�в�
            else        rtk->ssat[sat-1].resp[j/2]=v[nv];//α��в�
            
            /* variance���� */
            var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],j/2,j%2,opt)+
                    vart+SQR(C)*vari+var_rs[i];//varerr�������vart�����㷽�SQR(C)*vari����㷽�var_rs[i]����λ���Ӳ��õ������
            if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;
            
            trace(4,"%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n",str,sat,
                  j%2?"P":"L",j/2+1,v[nv],sqrt(var[nv]),azel[1+i*2]*R2D);
            
            /* reject satellite by pre-fit residuals ͨ��ǰ�в�ܾ����ǣ�������*/
            if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat,j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
                exc[i]=1; rtk->ssat[sat-1].rejc[j%2]++;
                continue;//�в���󣬾ܾ���������1��exc[i]=1��ʾ���ǲ�����
            }
            /* record large post-fit residuals ��¼�˲����в������*/
            if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {//�����ı������
                obsi[ne]=i; frqi[ne]=j; ve[ne]=v[nv]; ne++;//��¼��obsi[]��ʾ���ϸ����Ǳ�ţ�frqi[]��¼Ƶ�ʣ�ve[]�洢��вne��ʾ���
            }
            if (j%2==0) rtk->ssat[sat-1].vsat[j/2]=1;
            nv++;//nv����Ŀһֱ����������ֱ��n*nf-1��
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
        }//�ҵ��в��������ǡ�Ƶ�ʡ��в�ֵ
        sat=obs[maxobs].sat;
        trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
            post,str,sat,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[maxfrq%2]++; stat=0;
        ve[rej]=0;
    }
    /* constraint to local correction */
    nv+=const_corr(obs,n,exc,nav,x,pos,azel,rtk,v+nv,H+nv*rtk->nx,var+nv);
    
    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
        R[i+j*nv]=i==j?var[i]:0.0;//���
    }

	
	/*test R[]*/
	double c[100] = {0};
	int o, p;
	for (o = 0;o < 100;o++)
		for (p = 0;p < 100;p++){
			if (o == p) c[o] = R[o + o * nv];
		}

    return post?stat:nv;//post Ϊ�����stat���������nv
}
/* number of estimated states ����״̬��------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status ���½��״̬----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;
    
    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;
    
    if (rtk->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
    }
    rtk->sol.dtr[0]=rtk->x[IC(0,opt)];
    rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)];
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* test hold ambiguity ���Ա���ģ����-------------------------------------------------------*/
static int test_hold_amb(rtk_t *rtk)
{
    int i,j,stat=0;
    
    /* no fix-and-hold mode */
    if (rtk->opt.modear!=ARMODE_FIXHOLD) return 0;
    
    /* reset # of continuous fixed if new ambiguity introduced */
    for (i=0;i<MAXSAT;i++) {
        if (rtk->ssat[i].fix[0]!=2&&rtk->ssat[i].fix[1]!=2) continue;
        for (j=0;j<MAXSAT;j++) {
            if (rtk->ssat[j].fix[0]!=2&&rtk->ssat[j].fix[1]!=2) continue;
            if (!rtk->ambc[j].flags[i]||!rtk->ambc[i].flags[j]) stat=1;
            rtk->ambc[j].flags[i]=rtk->ambc[i].flags[j]=1;
        }
    }
    if (stat) {
        rtk->nfix=0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix>=rtk->opt.minfix;
}
/* precise point positioning ���ܵ��㶨λ-------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)//n -- the  satallite number
{
    const prcopt_t *opt=&rtk->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,dr[3]={0},std[3];//rs--���ǵ��ٶȺ�λ�ã�dts--�����Ӳ��������var--��λ���ı�����
    //v--����ֵ��ȥģ�ͻ�ֵ�òвH--��ƾ����ת�ã�R--��������Э���azel--��λ�Ǻ͸߶Ƚǣ�xp--δ֪������ֵ
	//Pp--���º��״̬��Э���dr--����ϫ����������std��׼ƫ�
	
	char str[32];
    int i,j,nv,info,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;//״̬ʱ����
    //svhΪ���ǵĽ���״��
   
	
	time2str(obs[0].time,str,2);//time2str����Ϊʱ���ʽ�����������ջ�ʱ��洢ep[]��str��
								//�����ʽΪyyyyMMddHHmmss.sssssssss   str���ڴ洢��׼ʱ��

    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n);//???
    
  
	
	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);//����������ڴ�ռ䣬mat���������ڴ�ռ�
    //rs--���ǵ��ٶȺ�λ�ã�dts--�����Ӳ��������var--������;n -- the  satallite number
    
	
	
	//����״̬�����������־��ʼ��
	for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;
    //i--number of satallite;j--number of freq; 
   
	
	
	/* temporal update of ekf states ��ǰ��չ�������˲�״̬���� */
    udstate_ppp(rtk,obs,n,nav);
    //״̬���£��������ջ�λ�ø���⼰Э������ջ��Ӳ��⼰Э��������㸡��⼰Э�������㣨���޵�������ʱ����ģ���ȵĸ���⼰Э���
    
	
	
	/* satellite positions and clocks����λ�����Ӳ� */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
	//ѡ����ʵ�������������ȡ���ǵ�λ�����ꡢ�ٶ��Լ������Ӳ����Ư�ƣ�varΪ����λ���Ӳ��
    
    if (rtk->opt.posopt[3]) {//����ѡ��Ϊ3
    /* exclude measurements of eclipsing satellite �ų��������ǵĲ���(block IIA) */
        testeclipse(obs,n,nav,rs);//n--����,rs--���ǵ�λ�ú��ٶ�
    }//�ų�������������Ϊ��


    /* earth tides correction ����ϫ����*/
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr==1?1:7,&nav->erp,opt->odisp[0],dr);
    }//opt->tidecorr��Ϊ��ʱ�����øú������г�ϫ��������������������dr[3]�С�
    nv=n*rtk->opt.nf*2+MAXSAT+3;//��λ��Ŀ+α����Ŀ+��������������ģ������Ŀ��+�Ӳ�+������+����㣻nvΪ����
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);//xpΪ����ֵ��RΪ��������
    v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);
    
    for (i=0;i<MAX_ITER;i++) {
        
        matcpy(xp,rtk->x,rtk->nx,1);
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
        
        /* prefit residuals  ǰ�в� */
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel))) {//��λ�ͱ���ʣ����xpΪδ֪��������ֵ����������װ��
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);//drΪ�����������v�в�
            break;
        }
        /* measurement update of ekf states ״̬��������*/
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {//��չ�������˲�
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }
        /* postfit residuals ��в�*/
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {//��λ�ͱ���ʣ��в�
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            stat=SOLQ_PPP;
            break;
        }
    }
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
    }
    if (stat==SOLQ_PPP) {
        
        /* ambiguity resolution in ppp */
        if (ppp_ar(rtk,obs,n,exc,nav,azel,xp,Pp)&&
            ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {
            
            matcpy(rtk->xa,xp,rtk->nx,1);
            matcpy(rtk->Pa,Pp,rtk->nx,rtk->nx);
            
            for (i=0;i<3;i++) std[i]=sqrt(Pp[i+i*rtk->nx]);
            if (norm(std,3)<MAX_STD_FIX) stat=SOLQ_FIX;//��Ϊģ���������ɵ��ж�����
        }
        else {
            rtk->nfix=0;//�̶���Ŀ
        }
        /* update solution status */
        update_stat(rtk,obs,n,stat);//���½���״̬
        
        /* hold fixed ambiguities */
        if (stat==SOLQ_FIX&&test_hold_amb(rtk)) {
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            trace(2,"%s hold ambiguity\n",str);
            rtk->nfix=0;
        }
    }
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
}
