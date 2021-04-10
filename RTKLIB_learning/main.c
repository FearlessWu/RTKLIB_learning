
#include "rtklib.h"

extern int showmsg(char *format, ...)
{
	return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {
}
extern void settime(gtime_t time) {
}
void main() {

	const prcopt_t prcopt_default = { /* defaults processing options */
		PMODE_SINGLE,0,2,SYS_GPS, /* mode,soltype,nf,navsys */
		15.0*D2R,{ { 0,0 } },         /* elmin,snrmask */
		EPHOPT_BRDC,1,1,1,                    /* sateph,modear,glomodear,bdsmodear */
		5,0,10,1,                   /* maxout,minlock,minfix,armaxiter */
		IONOOPT_IFLC,TROPOPT_SAAS,0,0,                    /* estion,esttrop,dynamics,tidecorr */
		1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
		0,0,                        /* rovpos,refpos */
		{ 100.0,100.0 },            /* eratio[] */
		{ 100.0,0.003,0.003,0.0,10.0 }, /* err[] */
		{ 30.0,0.03,0.3 },            /* std[] */
		{ 1E-4,1E-3,1E-4,1E-1,1E-2,0.0 }, /* prn[] */
		5E-12,                       /* sclkstab */
		{ 3.0,0.9999,0.25,0.1,0.05 },/* thresar */
		0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
		30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
		{ 0 },{ 0 },{ 0 },             /* baseline,ru,rb */
		{ "","" },                    /* anttype */
		{ { 0 } },{ { 0 } },{ 0 }     /* antdel,pcv,exsats */
	};

	const solopt_t solopt_default = { /* defaults solution output options */
		SOLF_XYZ,TIMES_GPST,1,3,    /* posf,times,timef,timeu */
		0,1,0,0,0,0,                /* degf,outhead,outopt,datum,height,geoid */
		0,0,0,                      /* solstatic,sstat,trace */
		{ 0.0,0.0 },                  /* nmeaintv */
		" ",""                      /* separator/program name */
	};
	filopt_t fopt = {
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
	};
	char *infile[] =
	{
		{ 0 },
		{ 0 },
		{ 0 },
		{ 0 },
	};
	int n;

	//infile[0] = "E:\\PPP\\BOS DATA\\tabl0010.16n";
	//infile[1] = "E:\\PPP\\BOS DATA\\tabl0010.16o";
	//infile[0] = "E:\\Examples\\2017244\\jfng2440.17o";
	//infile[0] = "..\\Examples\\2017244\\NNOR00AUS_R_20172440000_01D_30S_MO.rnx";
	//infile[1] = "..\\Examples\\2017244\\brdm2440.17p";
	//infile[2] = " "; 
	//infile[2] = "..\\Examples\\2017244\\wum19645.sp3";
	char *outfile = "..\\Examples\\2017244\\result\\NNOR.pos";
    infile[0] = "..\\Examples\\2017244\\MATE00ITA_R_20202220000_01D_30S_MO.20o";
    infile[1] = "..\\Examples\\2017244\\BRDC00IGS_R_20202220000_01D_MN.rnx";
    infile[2] = " ";
	n = 3;

	char *rov = " ", *base = " ";

	gtime_t ts = { 0 }, te = { 0 };
	double ti = 0.0, tu = 0.0;
	postpos(ts, te, ti, tu, &prcopt_default, &solopt_default, &fopt, infile, n, outfile, rov, base);
}