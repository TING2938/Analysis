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
int load_positionCenter(t_trxframe fr,int nm,int atoms,int nStart,double **pos,double Lbox[3],double *deloc, double totNCM);
int load_position(t_trxframe fr,int nm,int atoms,int nStart,double ***posmolecule,double Lbox[3]);
int load_velocity(t_trxframe fr,int nm,int atoms,int nStart,double ***velmolecule,double Lbox[3]);
int load_velocityCenter(t_trxframe fr,int nm,int atoms,int nStart,double **vel,double Lbox[3],double *deloc, double totNCM);
int vacf_compute(int nm, int index, double ***velcIndex, double *vacf);
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

double **** CreateMatrix_4d(int m,int n, int p, int q)
{
  int i;
  double ****m2;
  m2 = calloc((unsigned int) m, sizeof(double ***));
  for(i=0;i<m;i++)
    m2[i] = CreateMatrix_3d(n,p,q);
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
    "select which atom you want to examine with the -n argument.",
	"The result have been zoom in 10000, which should be scaled when final treatment"
  };
  
  /* default parameters */
  static int preframe  = 0;
  static int endframe  = 9999999;
  static int nStart1   = 0;            
  static int nm1       = 840;
  static int tN1       = 1;
  static int n_NCM     = 0; 
  static int dim       = 2;

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
	{ "-dim",FALSE, etINT, {&dim},
      "Dimension selected by user"
    }
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
  int        flags = TRX_READ_X | TRX_READ_V;   /* read position and velocity*/
  
  /* variables for this analysis */
  FILE    *inputPara,*outputfileR;
  int     i,j;  
  int     index,halfIndex,tmpNum,step,k,l,atoms1,atoms2,atoms_r,nType,m,n,sameMol,nij = 0,totN,numAt[6];
  double  atNCM[3][6][100],mol_NCM[3][6],dR[3],Lbox[3],L[3],dRc,dt,tmpt,label_r,label_icage,CR2;
  char    name0[60],temp[36];
  double  **posc1,**velc1,*Timeline,***pos1,***vel1,***velcIndex,vacf[99999],normalvacf[99999];
   
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
    for(i=0;i<numAt[j];i++){
      fscanf(inputPara,"%d%s%d%s%s%lf%lf%lf ",&atoms_r,temp,&atoms_r,temp,temp,&atNCM[0][j][i],&atNCM[1][j][i],&atNCM[2][j][i]);
      if(atNCM[0][j][i]>0) atNCM[0][j][i] = 1.0;
      else atNCM[0][j][i] = 0.0;
    }
  }
  tN1 --;  // make it as index in array

  atoms1 = numAt[tN1];
  // total mass/effective number/charge in a whole molecule (for solvent, the net charge is zero) 8/3/2010
  for(j=0;j<nType;j++){
    for(k=0;k<3;k++)
	{
      mol_NCM[k][j] = 0.0;
      for(i=0;i<numAt[j];i++)
	  {
		  mol_NCM[k][j] +=  atNCM[k][j][i];
	  }
      if(fabs(mol_NCM[k][j])<1e-6) 
	  {
		  mol_NCM[k][j] = 1.0;
	  }
    }
  }
  printf("\n# of atoms in reference type: %10d\n",atoms1);
  printf("Total:     %.1f %8.4f %8.4f \n",mol_NCM[0][tN1],mol_NCM[1][tN1],mol_NCM[2][tN1]);
  printf("           ---------------------\n");
  
  for(i=0;i<atoms1;i++)
    printf("           %.1f %8.4f %8.4f \n",atNCM[0][tN1][i],atNCM[1][tN1][i],atNCM[2][tN1][i]);

  /* GMX-467 support build-in function:
   * read_first_frame(oenv,&status,trx_file,&fr,flags);
   * read_next_frame(oenv,status,&fr);
   */
	double  ***dRcInd;
    printf("CreateMatrix.\n");
	posc1 = CreateMatrix(nm1,3);
	velc1 = CreateMatrix(nm1,3);
	Timeline     = CreateVector(100000);
	pos1 = CreateMatrix_3d(nm1, atoms1, 3);
	vel1 = CreateMatrix_3d(nm1, atoms1, 3);
	velcIndex = CreateMatrix_3d(nm1,3,99999);
	printf("CreateMatrix Done.\n");	
	read_first_frame(oenv,&status,trx_file,&fr,flags);
  
  
  for(step=0;step<preframe;step++)     
	  read_next_frame(oenv,status,&fr);
      
  tmpt = fr.time;
  /* =================== Main body of code ================== */
  //  initiate variables
  // type 1 & 2 are same, maybe atoms/groups different  
  index = 0;
  do{
    /*printf("Index = %d \n",index);*/	
    Lbox[0]   = fr.box[XX][XX]; Lbox[1]   = fr.box[YY][YY]; Lbox[2]   = fr.box[ZZ][ZZ];
	//load_position(fr, nm1, atoms1, nStart1, pos1, Lbox);
	//load_positionCenter(fr,nm1,atoms1,nStart1,posc1,Lbox,atNCM[n_NCM][tN1],mol_NCM[n_NCM][tN1]);	
	//load_velocity(fr, nm1, atoms1, nStart1, vel1, Lbox);
	load_velocityCenter(fr,nm1,atoms1,nStart1,velc1,Lbox,atNCM[n_NCM][tN1],mol_NCM[n_NCM][tN1]);
	for (i = 0; i < nm1; i++)
	{
		velcIndex[i][0][index] = velc1[i][0];
		velcIndex[i][1][index] = velc1[i][1];
		velcIndex[i][2][index] = velc1[i][2];
	}
	Timeline[index]=fr.time;
    dt  = tmpt;tmpt= fr.time;dt  = tmpt-dt;index ++;
  if(step==endframe) break;
  }while (read_next_frame(oenv,status,&fr));
  
  printf("Analysis done,4\n");
  
  //Calculate the VACF
  vacf_compute(nm1, index, velcIndex, vacf);
  
  halfIndex = index/2;
  for (i = 0; i < halfIndex; i++)
  {
	  vacf[i] = vacf[i]/halfIndex/nm1*10000;
	  //乘以10000是为了输出精度
  }
  
  for (i = 0; i < halfIndex; i++)
  {
	  normalvacf[i] = vacf[i] / vacf[0];
  }
  
  
  sprintf(name0, "VACF-mass-center-%d.dat",tN1+1);
  outputfileR = fopen(name0,"w");
  for (i = 0; i < halfIndex; i++)
  {
	  fprintf(outputfileR,"%8.2f  %8.4f  %8.4f\n",i*dt, vacf[i], normalvacf[i]);
  }
  fclose(outputfileR);
 
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

int load_velocity(t_trxframe fr,int nm,int atoms,int nStart,double ***velmolecule,double Lbox[3])
{
  int i,j,k,molN;
  double center[3], tmpVel[3],halfLbox[3];
  for(k=0;k<3;k++){halfLbox[k] = Lbox[k]/2;}
  for(i=0;i<nm;i++)
  {
    molN = nStart+i*atoms;
    for(j=0;j<atoms;j++){
		tmpVel[0] = fr.v[molN+j][XX];
		tmpVel[1] = fr.v[molN+j][YY];
		tmpVel[2] = fr.v[molN+j][ZZ];
		for(k=0;k<3;k++){velmolecule[i][j][k] = tmpVel[k];}
	}     
  }
  return 0;
}

int load_velocityCenter(t_trxframe fr,int nm,int atoms,int nStart,double **vel,double Lbox[3],double *deloc, double totNCM)
{
  int i,j,k,molN,m;
  double center[3], tmpVel, tmpCenter, halfLbox[3]; 
  for(i = 0 ; i < nm ; i++)
  {
    molN = nStart + i * atoms;  
    for(k = 0 ; k < 3 ; k++)
	{
      tmpCenter = 0.0;  // initiate for each molecule
      for(j = 0 ; j < atoms ; j++)
	  {
		tmpVel = fr.v[molN+j][k];	
		tmpCenter += tmpVel*deloc[j];
      }
      vel[i][k] = tmpCenter/totNCM;
    }
  }  
  return 0;	
}

int vacf_compute(int nm, int index, double ***velcIndex, double *vacf)
{
	int i,j,k,m,halfIndex;
	double tmp;
	halfIndex = index/2;
	for (i = 0; i < halfIndex; i++)
	{
		//数据初始化
		vacf[i] = 0.0;
	}
	for (i = 0; i < nm; i++)
	{
		for (j = 0; j < halfIndex; j++)
		{
			for (k = 0; k < halfIndex; k++)
			{
				tmp = 0;
				for (m = 0; m < 3; m++)
				{
					tmp = tmp + velcIndex[i][m][k+j]*velcIndex[i][m][k];
				}
				vacf[j] = vacf[j] + tmp ;
			}
		}
	}
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
