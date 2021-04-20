#include <string.h>
#include <stdlib.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/vec.h>
#include <gromacs/copyrite.h>
#include <gromacs/statutil.h>
#include <gromacs/tpxio.h>
#include <ctype.h>
#include <math.h>
#include <gromacs/sysstuff.h>
#include <gromacs/macros.h>
#include <gromacs/maths.h>
#include <gromacs/futil.h>
#include <gromacs/index.h>
#include <gromacs/xvgr.h>
#include <gromacs/gstat.h>
#include <gromacs/gmx_statistics.h>
#include <gromacs/pbc.h>
#include <gromacs/confio.h>
#include <gromacs/gmx_ana.h>
         
/* 
 *  get rdf of type1 to type2 based on molecules                                       2.10.2010
 *  add a parameter file for calculating the geometry center 
 *  calculate the m.s.d. of a type of molecules/atoms in a plane along the direction vertical to such plane 8.3.2010
 *  Coming from the code for MSD: calculate the RDF of a type of reference molecules/atoms with corresponding surrounding type
 *  in a plane along the direction vertical to such plane 8.12.2010
 *  note the bug:  if at1 and at2 in the same molecule, to get RDF of at1-at2, ideally, at2 should be at2 in other molecules
 *  i.e., at2 within the same molecule having at1 should not be accounted in.                               2.10.2011
 *  intermolecule or intramolecule                                                                          2.11.2011
 *  Support GMX-467 and have downward compatibility          4.21.2016 developed by Sheng Bi, mail address: chrishengbee@hust.edu.cn
 */
int msd_compute(int nm,int atoms, int halfchunkSize, int dim, int nbin, double *bin,
		double ***posCen, double **msd, double Lbox[3], double dbin, double lowPosition);
int load_alpha_epsilon(int atom1,int atom2,double **A,double **E,double *a1,double *a2, double *e1,double *e2);
int load_positionCenter(t_trxframe fr,int nm,int atoms,int nStart,double **pos,double Lbox[3],double *deloc, double totNCM);
int load_position(t_trxframe fr,int nm,int atoms,int nStart,double ***posmolecule,double Lbox[3]);
double periodicity(double dx,double box);
double * CreateVector(int cols);
double **CreateMatrix(int rows,int cols);
double *** CreateMatrix_3d(int m,int n, int p)
{
  int i;
  double ***m2;
  m2 = calloc((unsigned int) m, sizeof(double **));
  for(i=0;i<m;i++)
    m2[i] = CreateMatrix(n,p);
  return m2;
}

int main(int argc,char *argv[])
{
  /* need const char */
  const char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "when writing your own analysis tools. The advantage of ",
    "using gromacs for this is that you have access to all ",
    "information in the topology, and your program will be ",
    "able to handle all types of coordinates and trajectory ",
    "files supported by gromacs. Go ahead and try it! ",
    "This test version just writes the coordinates of an ",
    "arbitrary atom to standard out for each frame. You can ",
    "select which atom you want to examine with the -n argument."
  };
  
  /* default parameters */
  static int preframe  = 0;
  static int endframe  = 9999999;
  static int nStart1   = 0;        
  static int nStart2   = 10920;    
  static int nm1       = 840;
  static int nm2       = 840;
  static int tN1       = 1;
  static int tN2       = 2;
  static int n_NCM     = 0; 
  static int n_ncm     = 2;
  static int DIMN      = 2;
  static float dx      = 0.01;
  static float dy      = 0.01;
  static float dz      = 0.01;
  static float lowx    = 0.0;
  static float lowy    = 0.0;
  static float lowz    = 0.0;
  static float upx    = 5.0;
  static float upy    = 5.0;
  static float upz    = 10.0;
  static float cylCenter1 = 0.0;
  static float cylCenter2 = 0.0;
  static float dr     = 0.01;
  static float CR     = 0.5;
  static float Cr     = 0.0;

 printf("Test \n"); 
  
  /* GMX-467 support data structure */
  static const char *normtype[] = { NULL, "no", "x", "y", "z", NULL };
  static const char *axtitle[]  = { NULL, "no", "x", "y", "z", NULL };
  static gmx_bool    bTen       = FALSE;
  static gmx_bool    bMW        = TRUE;
  static gmx_bool    bRmCOMM    = FALSE;
  
  /* import parameters  */
   t_pargs pa[] = {
    { "-nstart1",FALSE, etINT, {&nStart1},
      "Start number of atoms in gro file"
    },
    { "-nm1",FALSE, etINT, {&nm1},
      "Number of molecules"
    },
    { "-tN1",FALSE, etINT, {&tN1},
      "the type of reference atom or molecule"
    },
    { "-NCM",FALSE, etINT, {&n_NCM},
      "center of molecule\n0: number; 1: charge; 2: mass"
    },
    { "-pre",FALSE, etINT, {&preframe},
      "Number of preframe"
    },
    { "-end",FALSE, etINT, {&endframe},
      "End number of trajectory set"
    },
    { "-ncm",FALSE, etINT, {&n_ncm},
      "# density\n0: number; 1: charge; 2: mass; when #=1,2 there is no center"
    },
	{ "-dim",FALSE, etINT, {&DIMN},
      "dim of bin sets\ndim>2 means getting D in 3-dimension; 0 means yz;1 means xz; 2 means xy;"
    },
	{ "-dx",FALSE, etREAL, {&dx},
      "the distance of dx(nm)"
    },
	{ "-dy",FALSE, etREAL, {&dy},
      "the distance of dy(nm)"
    },
	{ "-dz",FALSE, etREAL, {&dz},
      "the distance of dz(nm)"
    },
    { "-lowx",FALSE, etREAL, {&lowx},
      "the region of x(nm)"
    },
	{ "-lowy",FALSE, etREAL, {&lowy},
      "the region of y(nm)"
    },
	{ "-lowz",FALSE, etREAL, {&lowz},
      "the region of z(nm)"
    },
	{ "-upx",FALSE, etREAL, {&upx},
      "the region of x(nm)"
    },
	{ "-upy",FALSE, etREAL, {&upy},
      "the region of y(nm)"
    },
	{ "-upz",FALSE, etREAL, {&upz},
      "the region of z(nm)"
    },	
	{ "-c1",FALSE, etREAL, {&cylCenter1},
      "the center of the circle1 (nm)"
    },
	{ "-c2",FALSE, etREAL, {&cylCenter2},
      "the center of the circle2 (nm)"
    },
	{ "-dr",FALSE, etREAL, {&dr},
      "the dr bin of R (nm)"
    },
	{ "-CR",FALSE, etREAL, {&CR},
      "the up region of R (nm)"
    },
	{ "-Cr",FALSE, etREAL, {&Cr},
      "the low region of R (nm)"
    },	
	};

  /* GMX-467 support data structure */
  t_filenm           fnm[] = {
      { efTRX, NULL, NULL,  ffREAD },
      { efTPS, NULL, NULL,  ffREAD },
      { efNDX, NULL, NULL,  ffOPTRD },
      { efXVG, NULL, "msd", ffWRITE },
      { efXVG, "-mol", "diff_mol", ffOPTWR },
      { efPDB, "-pdb", "diff_mol", ffOPTWR }
  };
  
  /* GMX-467 support data structure */
  t_topology      top;
  int             ePBC;
  matrix          box;
  char            title[256];
  const char      *trx_file, *tps_file, *ndx_file, *msd_file, *mol_file, *pdb_file;
  rvec            *xdum;
  gmx_bool        bTop;
  int             axis, type;
  real            dim_factor;
  output_env_t    oenv;
  t_trxframe      fr;
  t_trxstatus     *status;
  int        flags = TRX_READ_X;   /* read only position */
  
  /* variables for this analysis */
  FILE    *inputPara,*outputfileR;
  int     i,j,k;  
  int     index,tmpNum,step,l,atoms1,atoms2,atoms_r,nType,m,n,sameMol,nij = 0,totN,numAt[6];
  double  atNCM[3][6][100],mol_NCM[3][6],dR[3],Lbox[3],L[3],dRc,dt,tmpt;
  char    name0[60],temp[36],MolType[60];
  double  **posc1,**posc2,***pos1,***pos2;
  


#define NFILE asize(fnm)

  parse_common_args(&argc, argv,
                    PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_BE_NICE,
                    NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
  trx_file = ftp2fn_null(efTRX, NFILE, fnm);
  tps_file = ftp2fn_null(efTPS, NFILE, fnm);
  bTop = read_tps_conf(tps_file, title, &top, &ePBC, &xdum, NULL, box, bMW || bRmCOMM);

  // get the atom information of geometry center
  // MSD_charge_para.dat comes from RDF_charge_para.dat
  if((inputPara = fopen("MSD_charge_para.dat","r")) == NULL){
    printf("\nNOTE: Parameter file cannot be loaded, mission terminated.\n");
    exit(0);
  }
  fscanf(inputPara,"%s %d \n",temp,&nType);
  if(tN1>nType){
    printf("\nNOTE: there is NO type %d in MSD_charge_para.dat.\n",tN1);
    exit(0);
  }
  for(j=0;j<nType;j++){
    fscanf(inputPara,"%s %d \n",temp,&numAt[j]);
	if (j == tN1-1) {strcpy(MolType,temp);}
    for(i=0;i<numAt[j];i++){
      fscanf(inputPara,"%d%s%d%s%s%lf%lf%lf ",&atoms_r,temp,&atoms_r,temp,temp,&atNCM[0][j][i],&atNCM[1][j][i],&atNCM[2][j][i]);
      if(atNCM[0][j][i]>0) atNCM[0][j][i] = 1.0;
      else atNCM[0][j][i] = 0.0;
    }
  }
  tN1 --; //tN2 --;  // make it as index in array

  atoms1 = numAt[tN1]; //atoms2 = numAt[tN2];
  // total mass/effective number/charge in a whole molecule (for solvent, the net charge is zero) 8/3/2010
  for(j=0;j<nType;j++){
    for(k=0;k<3;k++){
      mol_NCM[k][j] = 0.0;
      for(i=0;i<numAt[j];i++)
	mol_NCM[k][j] +=  atNCM[k][j][i];

      if(abs(mol_NCM[k][j])<1e-6) mol_NCM[k][j] = 1.0;
    }
  }
  printf("\n# of atoms in reference type: %10d\n",atoms1);
  printf("Total:     %.1f %8.4f %8.4f \n",mol_NCM[0][tN1],mol_NCM[1][tN1],mol_NCM[2][tN1]);
  printf("           ---------------------\n");
  
  for(i=0;i<atoms1;i++)
    printf("           %.1f %8.4f %8.4f \n",atNCM[0][tN1][i],atNCM[1][tN1][i],atNCM[2][tN1][i]);

/*
  printf("\n# of atoms in surrounding type: %10d\n",atoms2);
  printf("Total:     %.1f %8.4f %8.4f \n",mol_NCM[0][tN2],mol_NCM[1][tN2],mol_NCM[2][tN2]);
  printf("           ---------------------\n");
  
  for(i=0;i<atoms2;i++)
    printf("           %.1f %8.4f %8.4f \n",atNCM[0][tN2][i],atNCM[1][tN2][i],atNCM[2][tN2][i]);
*/
  /* GMX-467 support build-in function:
   * read_first_frame(oenv,&status,trx_file,&fr,flags);
   * read_next_frame(oenv,status,&fr);
   */
   
	double **NumDens1,**NumDens2,**CharDens1,**CharDens2;
	int IndexX,IndexY,IndexZ,nbinX,nbinY,nbinZ;
    printf("CreateMatrix.\n");
	posc1 = CreateMatrix(nm1,3);
	//posc2 = CreateMatrix(nm2,3);
	pos1  = CreateMatrix_3d(nm1,atoms1,3);
    //pos2  = CreateMatrix_3d(nm2,atoms2,3);
	printf("CreateMatrix Done.\n");
	


  /* =================== Main body of code ================== */
  //  initiate variables
  // type 1 & 2 are same, maybe atoms/groups different
  index = 0;
  nbinX=(upx-lowx)/dx;
  nbinY=(upy-lowy)/dy;
  nbinZ=(upz-lowz)/dz;
  double lowPos, upPos, lowD1, lowD2, upD1, upD2, dbin1, dbin2, tmpPos1, tmpPos2 ;
  int  nbin1, nbin2, me1, me2;
  double R, *NCMDens, *dV, tmpr, height;
  nbin1 = (CR-Cr)/dr;
  NCMDens = CreateVector(nbin1);
  dV = CreateVector(nbin1);
  
  if (DIMN == 0){height = upx - lowx; lowPos = lowx; upPos = upx;}
  if (DIMN == 1){height = upy - lowy; lowPos = lowy; upPos = upy;}  
  if (DIMN == 2){height = upz - lowz; lowPos = lowz; upPos = upz;}  
  for (i = 0; i < nbin1; i++)
  {
	  NCMDens[i] = 0.0;
	  tmpr = Cr + dr * (i+0.5);  
	  dV[i] = M_PI * (pow((i+1)*dr+Cr,2)-pow(i*dr+Cr,2)) * height;
  }
  
  read_first_frame(oenv,&status,trx_file,&fr,flags);  
  for(step=0;step<preframe;step++)     
	  read_next_frame(oenv,status,&fr);     
  tmpt = fr.time;
 
/*============================== number density ===================================*/  
if (n_ncm == 0)		//added 20190625
{
	do
	{
		Lbox[0]   = fr.box[XX][XX]; Lbox[1]   = fr.box[YY][YY]; Lbox[2]   = fr.box[ZZ][ZZ];
		load_position(fr,nm1,atoms1,nStart1,pos1,Lbox);
		load_positionCenter(fr,nm1,atoms1,nStart1,posc1,Lbox,atNCM[n_NCM][tN1],mol_NCM[n_NCM][tN1]);
		for (i = 0; i < nm1; i++)
		{
			if (DIMN == 0)
			{
				tmpPos1 = posc1[i][1] - cylCenter1;
				tmpPos2 = posc1[i][2] - cylCenter2;
			}
			if (DIMN == 1)
			{
				tmpPos1 = posc1[i][0] - cylCenter1;
				tmpPos2 = posc1[i][2] - cylCenter2;					
			}
			if (DIMN == 2)
			{
				tmpPos1 = posc1[i][0] - cylCenter1;
				tmpPos2 = posc1[i][1] - cylCenter2;
			}
			R = sqrt(tmpPos1*tmpPos1 + tmpPos2*tmpPos2);
			
			if (posc1[i][DIMN] >= lowPos && posc1[i][DIMN] <= upPos)
			{
				if (R >= Cr && R <= CR )
				{
					me1 = (R - Cr)/dr;
					NCMDens[me1]++;
				}								
			}
		}
		dt  = tmpt;tmpt= fr.time;dt  = tmpt-dt;index ++;	
		if(step==endframe) break;
	}while (read_next_frame(oenv,status,&fr));	
}
/*================================= charge & mass density ===========================*/
if (n_ncm == 1 || n_ncm == 2)		//added 20190625
{
	do
	{
		Lbox[0]   = fr.box[XX][XX]; Lbox[1]   = fr.box[YY][YY]; Lbox[2]   = fr.box[ZZ][ZZ];
		load_position(fr,nm1,atoms1,nStart1,pos1,Lbox);
		load_positionCenter(fr,nm1,atoms1,nStart1,posc1,Lbox,atNCM[n_NCM][tN1],mol_NCM[n_NCM][tN1]);
		for (i = 0; i < nm1; i++)
		{
			for (j = 0; j < atoms1; j++)
			{
				if (DIMN == 0)
				{
					tmpPos1 = pos1[i][j][1] - cylCenter1;
					tmpPos2 = pos1[i][j][2] - cylCenter2;
				}
				if (DIMN == 1)
				{
					tmpPos1 = pos1[i][j][0] - cylCenter1;
					tmpPos2 = pos1[i][j][2] - cylCenter2;
				}				
				if (DIMN == 2)
				{
					tmpPos1 = pos1[i][j][0] - cylCenter1;
					tmpPos2 = pos1[i][j][1] - cylCenter2;
				}
				R = sqrt(tmpPos1*tmpPos1 + tmpPos2*tmpPos2);				
				if (pos1[i][j][DIMN] >= lowPos && pos1[i][j][DIMN] <= upPos)
				{
					if (R >= Cr && R <= CR  )
					{
						me1 = (R - Cr)/dr;
						NCMDens[me1] = NCMDens[me1] + atNCM[n_ncm][tN1][j];
					}										
				}				
			}		
		}
		dt  = tmpt;tmpt= fr.time;dt  = tmpt-dt;index ++;	
		if(step==endframe) break;
	}while (read_next_frame(oenv,status,&fr));	
}

  sprintf(name0, "ncm_%d_Density_Circle-%.3f-%.3f_%s_bin1-%d_%1.2f-%1.2fnm-dim-%d.dat",n_ncm,cylCenter1,cylCenter2,MolType,nbin1,lowPos,upPos,DIMN);
  outputfileR = fopen(name0,"w");
  for(i=0;i<nbin1;i++)
  {
	  fprintf(outputfileR,"%10.3f %10.3f\n",(i+0.5)*dr+Cr,NCMDens[i]/dV[i]/index);
  }

   
  printf("\nTotal frame:    %10d\n",index);
  printf("Total time(ps): %10.3f\n",fr.time-tmpt);
  printf("Analysis done\n");
  thanx(stderr);

  return 0;
}

int msd_compute(int nm,int atoms, int halfchunkSize, int dim, int nbin, double *bin,
		double ***posCen, double **msd, double Lbox[3], double dbin, double lowPosition)
{
  int i,j,k,m,me,dimn;
  double r,r2;
  dimn = dim;
  if(dim>2) dim = 2;
  for(i=0;i<nm;i++){
    // for atom/molecule in the first time of a chunk, they will diffuse for a half chunk time
    me = (posCen[0][i][dim]-lowPosition)/dbin;
    if(me<nbin && me>-1){
      bin[me] += 1.0;
      
      for(k=0;k<halfchunkSize;k++)
	for(j=0;j<halfchunkSize;j++)
	  for(m=0;m<3;m++)
	    if(m!=dimn){
	      r  = periodicity(posCen[k+j][i][m]-posCen[j][i][m],Lbox[m]);
	      r2 = r*r;
	      msd[k][me] += r2;
	    }
    }
  }
  return 0;
}


int load_positionCenter(t_trxframe fr,int nm,int atoms,int nStart,double **pos,double Lbox[3],double *deloc, double totNCM)
{
  int i,j,k,molN,m;
  double center[3], tmpPos,tmpCenter,halfLbox[3];
  for(k=0;k<3;k++)
    halfLbox[k] = Lbox[k]/2;
  
  for(i=0;i<nm;i++){
    molN = nStart+i*atoms;
    
    for(k=0;k<3;k++){
      tmpCenter = 0.0;  // initiate for each molecule
      for(j=0;j<atoms;j++){
	tmpPos = fr.x[molN+j][k];
	// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
	while(tmpPos-fr.x[molN][k]> halfLbox[k])
	  tmpPos -= Lbox[k];
	while(tmpPos-fr.x[molN][k]<-halfLbox[k])
	  tmpPos += Lbox[k];
	
	tmpCenter += tmpPos*deloc[j];
      }
      pos[i][k] = tmpCenter/totNCM;
    }
  }
  
  return 0;
}

int load_position(t_trxframe fr,int nm,int atoms,int nStart,double ***posmolecule,double Lbox[3])
{
  int i,j,k,molN;
  double center[3], tmpPos[3],halfLbox[3];
  for(k=0;k<3;k++)
    halfLbox[k] = Lbox[k]/2;
  
  for(i=0;i<nm;i++){
    molN = nStart+i*atoms;
    for(j=0;j<atoms;j++)
      for(k=0;k<3;k++){
	tmpPos[k] = fr.x[molN+j][k];
	// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
	while(tmpPos[k]-fr.x[molN][k]> halfLbox[k])
	  tmpPos[k] -= Lbox[k];
	while(tmpPos[k]-fr.x[molN][k]<-halfLbox[k])
	  tmpPos[k] += Lbox[k];
	
	posmolecule[i][j][k] = tmpPos[k];
      }
  }

  return 0;
}
double periodicity(double dx,double box)
{
  while (dx > box/2.0)
    dx -=  box;
  while (dx < -box/2.0) 
    dx +=  box;

  return dx;
}

double * CreateVector(int rows)
{
   int	i;
   double  *m;

   m = calloc((unsigned int) rows,sizeof(double *));

   return m;
}
double ** CreateMatrix(int rows,int cols)
{
   int  i;
   double  **m;
 
   m = calloc((unsigned int) rows,sizeof(double *));
   for (i=0; i < rows; i++) {
      m[i] = calloc((unsigned int) cols,sizeof(double ));
   }
 
   return m;
}
