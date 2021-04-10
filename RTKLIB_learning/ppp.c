/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2016 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS 国际地球自转服务 tide model
*           -DOUTSTAT_AMB output ambiguity 模糊 parameters to solution status
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

#define MAX_ITER    8               /* max number of iterations 最大迭代次数 */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution 解决方案 */
#define MIN_NSAT_SOL 4              /* min satellite number for solution 解决方案的最少卫星数目 */
#define THRES_REJECT 4.0            /* reject threshold（阈值）of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance（变化） receiver position (m^2)开始变化接收机的位置 */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) 开始变化接收机时钟 */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2)开始变化ZTD的值 */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2)梯度 gradient（梯度）*/
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2)相位偏差 */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay 电离层延迟*/
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity（湿度） for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset 重置 ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular 测角的 velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

/* number and index of states 数量和状态指数  宏定义的函数*/
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)	//频率数量获取
#define NP(opt)     ((opt)->dynamics?9:3)						//动态模型
#define NC(opt)     (NSYS)										//导航系统数量获取
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))//对流层模型
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)		//电离层是否改正 改正模型
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))                           //9:3+NSYS
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)//3+1+(0/1/3)+i+1
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))//3+1+(0/1/3)+maxsat
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

/* standard deviation of state /浮点解 固定解的标准差 -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);//如果是固定状态那么。。。  nx--浮点数   SOLQ_FIX表示1为固定状态  sol.stat--求解类型  Pa固定解协方差
    return SQRT(rtk->P[i+i*rtk->nx]);// 浮点解协方差--P
}
/* write solution status for PPP /精密定位求解结果写出---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)//buff缓冲
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],*x;//
    int i,j,week;
    char id[32],*p=buff;
    
    if (!rtk->sol.stat) return 0;//sol_stat  求解类型为零时 
    
    trace(3,"pppoutstat:\n");// void trace(int leve,const char *format,....)?????
    
    tow=time2gpst(rtk->sol.time,&week);// void time2gpst(gtime_t ,char ,int)
    
    x=rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;//若固定，x=固定解，反之x=浮点解
    
    /* receiver position */
    p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2));
    
    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {            //opt.dynamics(0:none,1:velociy,2:accel)不等于0，--解算选项。动态方式
        ecef2pos(rtk->sol.rr,pos);      //sol.rr三维坐标与坐标改正数
        ecef2enu(pos,rtk->x+3,vel);    //rtk->x 浮点解
        ecef2enu(pos,rtk->x+6,acc);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                   "%.4f,%.5f,%.5f,%.5f\n",week,tow,rtk->sol.stat,vel[0],vel[1],
                   vel[2],acc[0],acc[1],acc[2],0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks  接收机时钟*/

    i=IC(0,&rtk->opt);
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,x[i]*1E9/CLIGHT,x[i+1]*1E9/CLIGHT,
               STD(rtk,i)*1E9/CLIGHT,STD(rtk,i+1)*1E9/CLIGHT);
    
    /* tropospheric parameters  对流层参数 */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {  //opt.tropopt--对流层选项，TROPOPT_EST=3 代表troposphere option: ZTD estimation
        i=IT(&rtk->opt);                                                  // TROPOPT_ESTG 4    troposphere option: ZTD+grad estimation 
        p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                   1,x[i],STD(rtk,i));
    }
    if (rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n",week,tow,
                   rtk->sol.stat,1,x[i+1],x[i+2],STD(rtk,i+1),STD(rtk,i+2));  /* double std [MAXSAT][3]-- fcb std-dev (cyc) */
    }
    /* ionosphere parameters 电离层参数 */
    if (rtk->opt.ionoopt==IONOOPT_EST) {   //#define IONOOPT_EST 4   ionosphere option: estimation 
        for (i=0;i<MAXSAT;i++) {        //MAXSAT--最大卫星数目
            ssat=rtk->ssat+i;            /* satellite status /卫星状态  包括方位角和高度*/
            if (!ssat->vs) continue;       /*ssat->vs=0 代表 valid satellite flag single /有效卫星标志？ */
            j=II(i+1,&rtk->opt);
            if (rtk->x[j]==0.0) continue;
            satno2id(i+1,id);
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,rtk->ssat[i].azel[0]*R2D,
                       rtk->ssat[i].azel[1]*R2D,x[j],STD(rtk,j));         /*FLOAT azel[2]--方位角和高度 azimuth/elevation (rad) */
        }
    }
#ifdef OUTSTAT_AMB       /* ambiguity parameters 模糊度参数*/
    
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
/* exclude meas of eclipsing satellite (block IIA) /排除消失/不显著/食卫星 ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");// void trace(int leve,const char *format,....)?????
    
    /* unit vector of sun direction (ecef)太阳方向的单位矢量 */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);//获得太阳与月亮的在地固坐标系下位置坐标
	//obs[0].time--接收机采样时间；erpv--地球自转参数；rsun--太阳地固坐标坐标
	/*void sunmoonpos(gtime_t tutc,const double *erpv, double *rsun,double *rmoon,double *gmst);   */

    normv3(rsun,esun);//计算太阳的方向余弦
	// int  normv3(const double *a, double *b)
    
    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;//antenna type
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;//距离不小于零
        
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;//从字符串str1中查找是否有字符串str2， 如果有，从str1中的str2位置起，返回str1的指针，如果没有，返回null
        
        /* sun-earth-satellite angle 太阳-地球-卫星角度*/
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);//cosa范围在-1到1之间，废话
        ang=acos(cosa);//获得卫星-太阳角度？？？，角度不大于pi/2；
        
        /* test eclipse  重叠 */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;//卫星的坐标置为零
    }
}
/* nominal yaw-angle /极小偏角处理 ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite /卫星偏角处理 return yaw ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model /卫星姿态模型  --------------------------------------------------*/
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
/* phase windup model 相位缠绕模型？？？-------------------------------------------------------*/
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
/* measurement error variance 测量方差------------------------------------------------*/
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
/* initialize state and covariance 初始状态和协方差-------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

/*x[]--浮点解；nx--浮点解个数，未知数个数，*/
/* geometry-free phase measurement 无几何相位测量-------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0) return 0.0;
    return lam[0]*obs->L[0]-lam[i]*obs->L[i];
}
/* Melbourne-Wubbena linear combination 线性组合 --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0||
        obs->P[0]==0.0||obs->P[i]==0.0) return 0.0;
    return lam[0]*lam[i]*(obs->L[0]-obs->L[i])/(lam[i]-lam[0])-
           (lam[i]*obs->P[0]+lam[0]*obs->P[i])/(lam[i]+lam[0]);
}
/* antenna corrected measurements  天线改正测量 --------------------------------------------*/
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
        
        /* antenna phase center and phase windup correction 天线相位中心和相位相缠改正*/
        L[i]=obs->L[i]*lam[i]-dants[i]-dantr[i]-phw*lam[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];
        
        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) 差分码偏差改正 */
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
    /* iono-free LC 无电离层模型*/
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
/* detect cycle slip by LLI 周跳检测 --------------------------------------------------*/
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
/* temporal update of position 暂存位置的更新-----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    int i;
    
    trace(3,"udpos_ppp:\n");
    
    /* fixed mode 固定模型*/
    if (rtk->opt.mode==PMODE_PPP_FIXED) {             // when the mode is ppp_fixedn
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);

		//生成以1E-8为协方差对角阵，矩阵大小为3*nx

        return;
    }//this mode is not been done yet.
	
	/* initialize position for first epoch 第一历元位置的初始化*/
    if (norm(rtk->x,3)<=0.0) {    //norm function :数组元素平方和开方，此处为距离小于0；
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);//receiver position vector is 60^2
    }


    /* static ppp mode静态ppp模型 */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {//模型为静态ppp的时候
        for (i=0;i<3;i++) {
            rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);//原来的协方差+噪声*历元时间差
        }
        return;
    }
	//当模型为静态ppp时，生成的协方差阵为P[i,i]=过程噪声^2*|delta t|


    /* kinmatic mode dynamics */
    for (i=0;i<3;i++) {
        initx(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
}
//receiver position vector is 60^2






/* temporal update of clock 时钟的时间更新--------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
    double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise)白噪声 每一历元钟初始化  */
    for (i=0;i<NSYS;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {  //当卫星星历为精密星历
            /* time of prec ephemeris is based gpst 精密星历是基于GPS时的*/
            /* negelect receiver inter-system bias 忽略接收机间的系统误差 */
            dtr=rtk->sol.dtr[3];//接收机钟差
			//dtr = rtk->sol.dtr[0];//接收机钟差
        }
        else {
            dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];//非精密星历
        }
        initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));//钟差浮点解初始化和协方差
//钟差浮点解为c*dtr,协方差为60的平方。
    }
}
/* temporal update of tropospheric parameters 对流层参数的时间更新--------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&rtk->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    
    if (rtk->x[i]==0.0) {//若对流层浮点解为零
        ecef2pos(rtk->sol.rr,pos);//地固坐标转到大地坐标，大地坐标格式为（latitude、longtitude、height）
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);//EGNOS模型对流层延迟改正
        initx(rtk,ztd,var,i);//对流层浮点解和协方差初始化，浮点解为对流层延迟改正数，协方差为对流层误差方差,矩阵大小为1*1
        
		if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(rtk,1E-6,VAR_GRA,j);
			//该种模式之下，再加两项，其浮点解为1E-6，协方差均为0.01^2,矩阵大小为2*2
        }
    }
    else {
        rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);//prn[2]--对流层噪声标准偏差值；得到协方差为：org+（std*t）^2
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);//prn[2]--对流层噪声标准偏差值；得到协方差为：org+（std*0.1*t）^2
            }
        }
    }
}
/* temporal update of ionospheric parameters 电离层的参数时间更新 ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double ion,sinel,pos[3],*azel;//lam载波波长、pos接收机位置、azel卫星天顶角
    char *p;
    int i,j,k,gap_resion=GAP_RESION;//重置电离层参数120
    
    trace(3,"udiono_ppp:\n");
    
    if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) {//strstr()函数是找出两个字符相同的位置，找到返回该位置指针，否则返回零。
        sscanf(p,"-GAP_RESION=%d",&gap_resion);
    }
    for (i=0;i<MAXSAT;i++) {
        j=II(i+1,&rtk->opt);
        if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) {//浮点解不等于零，且相位中断次数大于120，浮点解初始为0
            rtk->x[j]=0.0;
        }
    }
    for (i=0;i<n;i++) {
        j=II(obs[i].sat,&rtk->opt);
        if (rtk->x[j]==0.0) {
            k=satsys(obs[i].sat,NULL)==SYS_GAL?2:1;
            lam=nav->lam[obs[i].sat-1];//载波波长
            if (obs[i].P[0]==0.0||obs[i].P[k]==0.0||lam[0]==0.0||lam[k]==0.0) {//P为 伪距观测数据
                continue;
            }
            ion=(obs[i].P[0]-obs[i].P[k])/(1.0-SQR(lam[k]/lam[0]));
            ecef2pos(rtk->sol.rr,pos);//地固坐标转大地坐标
            azel=rtk->ssat[obs[i].sat-1].azel;//方位角数组
            ion/=ionmapf(pos,azel);//电离层穿透点（IPP）位置和倾斜因子
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
        initx(rtk,1E-6,VAR_DCB,i);//浮点解为1E-6，协方差为30平方。
    }
}
/* temporal update of phase biases 相位偏差-------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double ion,dantr[NFREQ]={0},dants[NFREQ]={0};
    int i,j,k,l,f,sat,slip[MAXOBS]={0},clk_jump=0;
    
    trace(3,"udbias  : n=%d\n",n);
    
    /* handle day-boundary clock jump//处理日边界时钟跳转 */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;//初始化周跳标识
    }
    /* detect cycle slip by LLI 周跳探测 */
    detslp_ll(rtk,obs,n);
    
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk,obs,n,nav);
    
    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(rtk,obs,n,nav);
//周跳检测 模型，如果存在周跳：rtk->ssat[obs[i].sat-1].slip[j]|=1
    
    ecef2pos(rtk->sol.rr,pos);
    
    for (f=0;f<NF(&rtk->opt);f++) {
        
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));//如果观测数据的中断计数器失效，重置浮点解与协方差阵
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);
            corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,dantr,dants,
                      0.0,L,P,&Lc,&Pc);//天线改正测量（改正了L、P的观测值）
            
            bias[i]=0.0;
            
            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc-Pc;
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];//无电离层组合的选项
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
/* temporal update of states 状态的更新--------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(3,"udstate_ppp: n=%d\n",n);
    
    /* temporal update of position位置更新 */
    udpos_ppp(rtk);
	
/*//位置更新包括接收机三维坐标更新与浮点解的协方差阵的更新
when the mode is ppp_fixed，生成以1E-8为协方差对角阵，矩阵大小为3*nx,第一历元位置的初始化,
 receiver position vector is 60^2;
模型为静态ppp的时,生成的协方差阵为P[i,i]=过程噪声^2*|delta t|;
动态时，receiver position vector is 60^2 */
    

    /* temporal update of clock 钟差更新*/
    udclk_ppp(rtk);
    /*钟差更新：钟差浮点解为c*dtr,协方差为60的平方，矩阵大小为[NSYS*NSYS]*/

    /* temporal update of tropospheric parameters 对流层参数更新 */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
	/*对流层更新：
	当对流层选项为Estimate ZTD或者Estimate ZTD+Grad
	浮点解为对流层延迟改正数ztd（EGNOS模型）；
	原来浮点解为0时，即x[i]=0:
	                         协方差为对流层误差方差var,矩阵大小为1*1
	 对流层选项为TROPOPT_ESTG时，再加两项，其浮点解为1E-6，协方差均为0.01^2,矩阵大小为2*2
	else：
	     得到协方差为：org+（std*t）^2
		 对流层选项为TROPOPT_ESTG时，协方差阵再加两项，得到协方差为：org+（std*0.0*t）^2
	*/



    /* temporal update of ionospheric parameters 电离层参数更新 */
    if (rtk->opt.ionoopt==IONOOPT_EST) {
        udiono_ppp(rtk,obs,n,nav);
    }


    /* temporal update of L5-receiver-dcb parameters码间偏差参数 */
    if (rtk->opt.nf>=3) {  //当频率数大于3的时
        uddcb_ppp(rtk);//浮点解为1E-6，协方差为30平方
    }


    /* temporal update of phase-bias //相位偏差更新*/
    udbias_ppp(rtk,obs,n,nav);//模糊度的浮点解=bais[i]以及协方差60square

}

/* satellite antenna phase center variation 卫星天线相位中心变量----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;
    
    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;//eu为卫星方向余弦；ez为接收机以地心为原点的方向余弦
    
    cosa=dot(eu,ez,3);//
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);//角度
    
    antmodel_s(pcv,nadir,dant);//天线参数变量dant
}
/* precise tropospheric model 精确对流层模型------------------------------------------------*/
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
    
    if (opt->tropopt==TROPOPT_SAAS) {//萨斯塔莫宁模型
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);//dtrp为其对流层延迟
        *var=SQR(ERR_SAAS);//var为方差
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {//SBAS模型
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {//ZTD/ZTD+GRAD估计
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);//精密对流层模型
        return 1;
    }
    if (opt->tropopt==TROPOPT_ZTD) {
        if (pppcorr_trop(&nav->pppcorr,time,pos,trp,std)) {
            *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);//输入改正参数
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
/* constraint to local correction 局部修正约束--------------------------------------------*/
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
        pppcorr_trop(&nav->pppcorr,time,pos,trop,std_trop)) {//对流层改正
        
        for (i=0;i<(rtk->opt.tropopt==TROPOPT_EST?1:3);i++) {
            if (std_trop[i]==0.0) continue;
            j=IT(&rtk->opt)+i;
            v[nv]=trop[i]-x[j];//对流层改正数
            for (k=0;k<rtk->nx;k++) H[k+nv*rtk->nx]=k==j?1.0:0.0;
            var[nv++]=SQR(std_trop[i]);//对流层方差
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
/* phase and code residuals 相位和编码剩余误差--------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel)//rs--卫星位置与速度；dts--卫星钟差；var_rs--卫星位置钟差误差方差；svh--卫星健康状况；dr--潮汐坐标改正数
	//exc--卫星是否排除；x--未知数近似值；v--？？；H--设计矩阵的转置；R--测量噪声
{
    const double *lam;//波长
    prcopt_t *opt=&rtk->opt;
    double y,r,cdtr,bias,C,rr[3],pos[3],e[3],dtdx[3],L[NFREQ],P[NFREQ],Lc,Pc;
    double var[MAXOBS*2],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,dcb;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    double ve[MAXOBS*2*NFREQ]={0},vmax=0;
    char str[32];
    int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
    int i,j,k,sat,sys,nv=0,nx=rtk->nx,stat=1;
    
    time2str(obs[0].time,str,2);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;//对每个卫星进行有效卫星标志进行初始化，初始化值为0，频率数最多3个
    
    for (i=0;i<3;i++) rr[i]=x[i]+dr[i];//dr为潮汐改正参数
    ecef2pos(rr,pos);//接收机地固坐标rr转到大地坐标系pos
    
    for (i=0;i<n&&i<MAXOBS;i++) {//对每个卫星进行操作
        sat=obs[i].sat;//satellite number
        lam=nav->lam[sat-1];//载波波长
        if (lam[j/2]==0.0||lam[0]==0.0) continue;//j=2；当波长为零跳出本次循环
        
        if ((r=geodist(rs+i*6,rr,e))<=0.0||//rs+i*6表示数组第i*6+1个元素；此函数获得接收机到卫星的方向余弦，r为加上地球自转参数的后的距离
            satazel(pos,e,azel+i*2)<opt->elmin) {//获得卫星的方位角和高度角azel；e，地固坐标下的接收机向量
            exc[i]=1;//排除卫星
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||//1.转换卫星号获取卫星系统 2.有效卫星标志
            satexclude(obs[i].sat,svh[i],opt)||exc[i]) {//3.排除高度角不符合的卫星
            exc[i]=1;
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&vart)||//获得对流层改正数和协方差
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
            continue;
        }
        /* satellite and receiver antenna model */
        if (opt->posopt[0]) satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);//dants卫星天线参数变量
        antmodel(opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr);//dantr接收机天线改正量
        
        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {//相位缠绕改正
            continue;
        }
        /* corrected phase and code measurements 相位与码测量改正*/
        corr_meas(obs+i,nav,azel+i*2,&rtk->opt,dantr,dants,
                  rtk->ssat[sat-1].phw,L,P,&Lc,&Pc);
        
        /* stack phase and code residuals {L1,P1,L2,P2,...}堆砌相位与编码残差 */
        for (j=0;j<2*NF(opt);j++) {
            
            dcb=bias=0.0;
            
            if (opt->ionoopt==IONOOPT_IFLC) {//ionosphere option: L1/L2 or L1/L5 iono-free LC 
                if ((y=j%2==0?Lc:Pc)==0.0) continue;//j=0,2为载波相位
            }
            else {
                if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;
                
                /* receiver DCB correction for P2 */
                if (j/2==1) dcb=-nav->rbias[0][sys==SYS_GLO?1:0][0];//rbias[64][2][3];接收机码偏差
            }
            C=SQR(lam[j/2]/lam[0])*ionmapf(pos,azel+i*2)*(j%2==0?-1.0:1.0);//pos接收机大地坐标
            for (k=0;k<nx;k++) H[k+nx*nv]=k<3?-e[k]:0.0;//对前三个的系数矩阵，为坐标向量的得相反数
            
            /* receiver clock */
            k=sys==SYS_GLO?1:(sys==SYS_GAL?2:(sys==SYS_CMP?3:0));
            cdtr=x[IC(k,opt)];//接收机钟差粗值
            H[IC(k,opt)+nx*nv]=1.0;//钟差系数矩阵
            
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
                H[IB(sat,j/2,opt)+nx*nv]=1.0;//IB的位置在变，保证了模糊度在矩阵的位置
            }
            /* residual */
            v[nv]=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias);//相位，伪距的残差
            
            if (j%2==0) rtk->ssat[sat-1].resc[j/2]=v[nv];//相位残差
            else        rtk->ssat[sat-1].resp[j/2]=v[nv];//伪距残差
            
            /* variance方差 */
            var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],j/2,j%2,opt)+
                    vart+SQR(C)*vari+var_rs[i];//varerr测量方差；vart对流层方差；SQR(C)*vari电离层方差；var_rs[i]卫星位置钟差方差；得到方差和
            if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;
            
            trace(4,"%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n",str,sat,
                  j%2?"P":"L",j/2+1,v[nv],sqrt(var[nv]),azel[1+i*2]*R2D);
            
            /* reject satellite by pre-fit residuals 通过前残差拒绝卫星？？？？*/
            if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat,j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
                exc[i]=1; rtk->ssat[sat-1].rejc[j%2]++;
                continue;//残差过大，拒绝计数器加1，exc[i]=1表示卫星不包括
            }
            /* record large post-fit residuals 记录滤波后大残差？？？？*/
            if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {//大于四倍中误差
                obsi[ne]=i; frqi[ne]=j; ve[ne]=v[nv]; ne++;//记录，obsi[]表示不合格卫星编号；frqi[]记录频率；ve[]存储大残差，ne表示序号
            }
            if (j%2==0) rtk->ssat[sat-1].vsat[j/2]=1;
            nv++;//nv的数目一直会往下增加直到n*nf-1；
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
        }//找到残差最大的卫星、频率、残差值
        sat=obs[maxobs].sat;
        trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
            post,str,sat,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[maxfrq%2]++; stat=0;
        ve[rej]=0;
    }
    /* constraint to local correction */
    nv+=const_corr(obs,n,exc,nav,x,pos,azel,rtk,v+nv,H+nv*rtk->nx,var+nv);
    
    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
        R[i+j*nv]=i==j?var[i]:0.0;//误差
    }

	
	/*test R[]*/
	double c[100] = {0};
	int o, p;
	for (o = 0;o < 100;o++)
		for (p = 0;p < 100;p++){
			if (o == p) c[o] = R[o + o * nv];
		}

    return post?stat:nv;//post 为真输出stat，否则输出nv
}
/* number of estimated states 估计状态数------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status 更新解决状态----------------------------------------------------*/
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
/* test hold ambiguity 测试保留模糊度-------------------------------------------------------*/
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
/* precise point positioning 精密单点定位-------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)//n -- the  satallite number
{
    const prcopt_t *opt=&rtk->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,dr[3]={0},std[3];//rs--卫星的速度和位置；dts--卫星钟差与改正；var--相位中心变量；
    //v--测量值减去模型化值得残差；H--设计矩阵的转置；R--测量误差的协方差；azel--方位角和高度角；xp--未知数近似值
	//Pp--更新后的状态的协方差；dr--地球潮汐改正参数，std标准偏差；
	
	char str[32];
    int i,j,nv,info,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;//状态时单点
    //svh为卫星的健康状况
   
	
	time2str(obs[0].time,str,2);//time2str函数为时间格式函数，将接收机时间存储ep[]在str中
								//输出格式为yyyyMMddHHmmss.sssssssss   str用于存储标准时间

    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n);//???
    
  
	
	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);//给矩阵分配内存空间，mat函数分配内存空间
    //rs--卫星的速度和位置；dts--卫星钟差与改正；var--？？？;n -- the  satallite number
    
	
	
	//卫星状态的整周数解标志初始化
	for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;
    //i--number of satallite;j--number of freq; 
   
	
	
	/* temporal update of ekf states 当前扩展卡尔曼滤波状态更新 */
    udstate_ppp(rtk,obs,n,nav);
    //状态更新：包括接收机位置浮点解及协方差、接收机钟差浮点解及协方差、对流层浮点解及协方差、电离层（非无电离层组合时）、模糊度的浮点解及协方差。
    
	
	
	/* satellite positions and clocks卫星位置与钟差 */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
	//选择合适的卫星星历，获取卫星的位置坐标、速度以及卫星钟差和钟漂移，var为卫星位置钟差方差
    
    if (rtk->opt.posopt[3]) {//解算选项为3
    /* exclude measurements of eclipsing satellite 排除？？卫星的测量(block IIA) */
        testeclipse(obs,n,nav,rs);//n--卫星,rs--卫星的位置和速度
    }//排除的卫星坐标置为零


    /* earth tides correction 地球潮汐改正*/
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr==1?1:7,&nav->erp,opt->odisp[0],dr);
    }//opt->tidecorr不为零时，调用该函数进行潮汐改正。坐标改正数存放在dr[3]中。
    nv=n*rtk->opt.nf*2+MAXSAT+3;//相位数目+伪距数目+最大卫星数（最大模糊度数目）+钟差+对流层+电离层；nv为列数
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);//xp为测量值，R为测量方差
    v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);
    
    for (i=0;i<MAX_ITER;i++) {
        
        matcpy(xp,rtk->x,rtk->nx,1);
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
        
        /* prefit residuals  前残差 */
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel))) {//相位和编码剩余误差，xp为未知参数近似值，将矩阵组装好
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);//dr为坐标改正数，v残差
            break;
        }
        /* measurement update of ekf states 状态测量更新*/
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {//扩展卡尔曼滤波
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }
        /* postfit residuals 后残差*/
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {//相位和编码剩余残差
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
            if (norm(std,3)<MAX_STD_FIX) stat=SOLQ_FIX;//认为模糊度求解完成的判断条件
        }
        else {
            rtk->nfix=0;//固定数目
        }
        /* update solution status */
        update_stat(rtk,obs,n,stat);//更新解算状态
        
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
