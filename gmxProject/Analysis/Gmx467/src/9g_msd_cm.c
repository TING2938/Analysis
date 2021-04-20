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
		double ***posCen, double **msd, double Lbox[3], double dbin, double lowPosition, int region);
int load_positionCenter(t_trxframe fr,int nm,int atoms,int nStart,double **pos,double Lbox[3],double *deloc, double totNCM);
int load_position(t_trxframe fr,int nm,int atoms,int nStart,double ***posmolecule,double Lbox[3]);
double periodicity(double dx,double box);
double * CreateVector(int cols);
double **CreateMatrix(int rows,int cols);
double*** CreateMatrix_3d(int m, int n, int p)
{
	int i;
	double*** m2;
	m2 = (double***) calloc((unsigned int)m, sizeof(double**));
	for (i = 0; i < m; i++)
		m2[i] = CreateMatrix(n, p);
	return m2;
}

int main(int argc, char* argv[])
{
	/* need const char */
	const char* desc[] = {
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
	static int preframe = 0;
	static int endframe = 9999999;
	static int nStart1 = 0;
	static int nm1 = 840;
	static int tN1 = 1;
	static int n_NCM = 0;
	static int dim = 2;
	static float lowPos = 0.0;
	static float upPos = 10.0;
	static int nbin = 1;
	static int region = 1;

	printf("Test \n");

	/* GMX-467 support data structure */
	static const char* normtype[] = { NULL, "no", "x", "y", "z", NULL };
	static const char* axtitle[] = { NULL, "no", "x", "y", "z", NULL };
	static gmx_bool    bTen = FALSE;
	static gmx_bool    bMW = TRUE;
	static gmx_bool    bRmCOMM = FALSE;

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
	 },
	 { "-nbin",FALSE, etINT, {&nbin},
	   "nbin selected by user"
	 },

	 { "-lowPos",FALSE, etREAL, {&lowPos},
	   "the region of lowPos(nm)"
	 },
	 { "-upPos",FALSE, etREAL, {&upPos},
	   "the region of upPos(nm)"
	 },
	 { "-region",FALSE, etINT, {&region},
	   "the region selected by user; 1 means X; 2 means Y; 3 means Z; 4 means XY; 5 means XYZ"
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
	const char* trx_file, * tps_file, * ndx_file, * msd_file, * mol_file, * pdb_file;
	rvec* xdum;
	gmx_bool        bTop;
	int             axis, type;
	real            dim_factor;
	output_env_t    oenv;
	t_trxframe      fr;
	t_trxstatus* status;
	int        flags = TRX_READ_X;   /* read only position */

	/* variables for this analysis */
	FILE* inputPara, * outputfileR;
	int     i, j, k;
	int     index, tmpNum, step, l, atoms1, atoms_r, nType, m, n, sameMol, nij = 0, totN, numAt[6];
	double  atNCM[3][6][100], mol_NCM[3][6], dR[3], Lbox[3], L[3], dRc, dt, tmpt;
	char    name0[60], temp[36];
	double*** posc1;



#define NFILE asize(fnm)

	parse_common_args(&argc, argv,
		PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_BE_NICE,
		NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	trx_file = ftp2fn_null(efTRX, NFILE, fnm);
	tps_file = ftp2fn_null(efTPS, NFILE, fnm);
	bTop = read_tps_conf(tps_file, title, &top, &ePBC, &xdum, NULL, box, bMW || bRmCOMM);

	// get the atom information of geometry center
	// MSD_charge_para.dat comes from RDF_charge_para.dat
	if ((inputPara = fopen("MSD_charge_para.dat", "r")) == NULL) {
		printf("\nNOTE: Parameter file cannot be loaded, mission terminated.\n");
		exit(0);
	}
	fscanf(inputPara, "%s %d \n", temp, &nType);
	if (tN1 > nType) {
		printf("\nNOTE: there is NO type %d in MSD_charge_para.dat.\n", tN1);
		exit(0);
	}
	for (j = 0; j < nType; j++) {
		fscanf(inputPara, "%s %d \n", temp, &numAt[j]);
		for (i = 0; i < numAt[j]; i++) {
			fscanf(inputPara, "%d%s%d%s%s%lf%lf%lf ", &atoms_r, temp, &atoms_r, temp, temp, &atNCM[0][j][i], &atNCM[1][j][i], &atNCM[2][j][i]);
			if (atNCM[0][j][i] > 0) atNCM[0][j][i] = 1.0;
			else atNCM[0][j][i] = 0.0;
		}
	}
	tN1--;   // make it as index in array

	atoms1 = numAt[tN1];
	// total mass/effective number/charge in a whole molecule (for solvent, the net charge is zero) 8/3/2010
	for (j = 0; j < nType; j++) {
		for (k = 0; k < 3; k++) {
			mol_NCM[k][j] = 0.0;
			for (i = 0; i < numAt[j]; i++)
				mol_NCM[k][j] += atNCM[k][j][i];

			if (fabs(mol_NCM[k][j]) < 1e-6) mol_NCM[k][j] = 1.0;
		}
	}
	printf("\n# of atoms in reference type: %10d\n", atoms1);
	printf("Total:     %.1f %8.4f %8.4f \n", mol_NCM[0][tN1], mol_NCM[1][tN1], mol_NCM[2][tN1]);
	printf("           ---------------------\n");
	for (i = 0; i < atoms1; i++)
		printf("           %.1f %8.4f %8.4f \n", atNCM[0][tN1][i], atNCM[1][tN1][i], atNCM[2][tN1][i]);

	/* GMX-467 support build-in function:
	 * read_first_frame(oenv,&status,trx_file,&fr,flags);
	 * read_next_frame(oenv,status,&fr);
	 */
	int halfIndex, DIMN;
	double dbin, dx, * numbin, ** msdCenter, ** posctmp;
	printf("CreateMatrix.\n");
	posc1 = CreateMatrix_3d(49999, nm1, 3);
	printf("CreateMatrix Done.\n");
	dbin = (upPos - lowPos) / nbin;
	DIMN = dim;
	numbin = CreateVector(nbin);
	msdCenter = CreateMatrix(49999, nbin);
	posctmp = CreateMatrix(nm1, 3);

	read_first_frame(oenv, &status, trx_file, &fr, flags);

	for (step = 0; step < preframe; step++)
		read_next_frame(oenv, status, &fr);

	tmpt = fr.time;
	/* =================== Main body of code ================== */
	//  initiate variables
	// type 1 & 2 are same, maybe atoms/groups different
	index = 0;
	do {
		/*printf("Index = %d \n",index);*/
		// printf("\n");
		Lbox[0] = fr.box[XX][XX]; 
		Lbox[1] = fr.box[YY][YY]; 
		Lbox[2] = fr.box[ZZ][ZZ];
		load_positionCenter(fr, nm1, atoms1, nStart1, posctmp, Lbox, atNCM[n_NCM][tN1], mol_NCM[n_NCM][tN1]);
		for (i = 0; i < nm1; i++)
		{
			for (j = 0; j < 3; j++)
			{
				posc1[index][i][j] = posctmp[i][j];
			}
			// 	printf("i = %d pos = %f %f %f\n",i,posc1[index][i][0],posc1[index][i][1],posc1[index][i][2]);
		}

		dt = tmpt; 
		tmpt = fr.time; 
		dt = tmpt - dt; 
		index++;
		if (step == endframe) 
			break;
	} while (read_next_frame(oenv, status, &fr));
	printf("\ntime step (ps):  %10.3f\n", dt);
	printf("\nTotal frame:    %10d\n", index);

	for (i = 1; i < index; i++)
	{
		for (k = 0; k < nm1; k++)
		{
			for (m = 0; m < 3; m++)
			{
				//adjust the periodicity
				while ((posc1[i][k][m] - posc1[i - 1][k][m]) > Lbox[m] / 2) 
				{ 
					posc1[i][k][m] = posc1[i][k][m] - Lbox[m]; 
				}
				while ((posc1[i][k][m] - posc1[i - 1][k][m]) <= -Lbox[m] / 2) 
				{ 
					posc1[i][k][m] = posc1[i][k][m] + Lbox[m]; 
				}
			}
			//  printf("i = %d pos = %f %f %f\n",i,posc1[i][k][0],posc1[i][k][1],posc1[i][k][2]);
		}
	}
	halfIndex = index / 2;
	msd_compute(nm1, atoms1, halfIndex, DIMN, nbin, numbin, posc1, msdCenter, Lbox, dbin, lowPos, region);
	if (region == 1)
	{
		sprintf(name0, "MSD-center_bin%d_%dtype_%1.1f_%1.1f_X.dat", nbin, tN1 + 1, lowPos, upPos);
	}
	if (region == 2)
	{
		sprintf(name0, "MSD-center_bin%d_%dtype_%1.1f_%1.1f_Y.dat", nbin, tN1 + 1, lowPos, upPos);
	}
	if (region == 3)
	{
		sprintf(name0, "MSD-center_bin%d_%dtype_%1.1f_%1.1f_Z.dat", nbin, tN1 + 1, lowPos, upPos);
	}
	if (region == 4)
	{
		sprintf(name0, "MSD-center_bin%d_%dtype_%1.1f_%1.1f_XY.dat", nbin, tN1 + 1, lowPos, upPos);
	}
	if (region == 5)
	{
		sprintf(name0, "MSD-center_bin%d_%dtype_%1.1f_%1.1f_XYZ.dat", nbin, tN1 + 1, lowPos, upPos);
	}


	printf("The output name is: %s! for m.s.d based on molecule center!!\n", name0);
	printf("nbin_AtomType_chunkTime_bigenTime!!\n");
	printf("Analysis done\n");

	outputfileR = fopen(name0, "w");
	for (i = 0; i < halfIndex; i++)
	{
		fprintf(outputfileR, "%8.3f  ", i * dt);
		for (j = 0; j < nbin; j++)
		{
			//printf("numbin %f msd %f halfIndex %d\n",numbin[j],msdCenter[i][j],halfIndex);	  
			fprintf(outputfileR, "%10.6f  ", msdCenter[i][j]);
		}
		fprintf(outputfileR, "\n");
	}

	printf("\nTotal frame:    %10d\n", index);
	printf("Total time(ps): %10.3f\n", fr.time - tmpt);
	printf("Analysis done\n");
	thanx(stderr);

	return 0;
}

int msd_compute(int nm, int atoms, int halfIndex, int dim, int nbin, double* bin,
	double*** posCen, double** msd, double Lbox[3], double dbin, double lowPosition, int region)
{
	int i, j, k, m, l, me, dimn;
	double r[3] = { 0.0,0.0,0.0 }, r2 = 0.0, tmp = 0.0;
	dimn = dim;
	if (dim > 2) dim = 2;

	int Dim_factor;
	for (i = 0; i < nbin; i++) 
	{ 
		bin[i] = 0; 
	}
	for (i = 0; i < halfIndex; i++)
	{
		for (j = 0; j < nbin; j++)
		{
			msd[i][j] = 0.0;
		}
	}
	if (region == 1 || region == 2 || region == 3)
	{
		//calculate the X or Y or Z msd
		region = region - 1;
		Dim_factor = 2;
		for (i = 0; i < nm; i++)
		{
			for (k = 0; k < halfIndex; k++)
			{
				tmp = posCen[k][i][dim];    // calculate the bin by dim;
				while (tmp > Lbox[dim]) 
				{ 
					tmp = tmp - Lbox[dim]; 
				}
				while (tmp < 0) 
				{ 
					tmp = tmp + Lbox[dim]; 
				}
				me = (tmp - lowPosition) / dbin;
				//me=(posCen[k][i][dim]-lowPosition)/dbin;
				if (me < nbin && me > -1)
				{
					bin[me] = bin[me] + 1.0;
					for (j = 0; j < halfIndex; j++)
					{
						// r=periodicity(posCen[k+j][i][region]-posCen[k][i][region],Lbox[region]);
						l = k + j;
						r[region] = posCen[l][i][region] - posCen[k][i][region];
						r2 = r[region] * r[region];
						msd[j][me] = msd[j][me] + r2;
						// 	printf(" l %d k %d posl %f posk %f msd %f\n",l, k,posCen[l][i][region],posCen[k][i][region], msd[j][me]);
					}
				}
			}
		}
	}
	if (region == 4)
	{
		//calculate the XY msd
		Dim_factor = 4;
		for (i = 0; i < nm; i++)
		{
			for (k = 0; k < halfIndex; k++)
			{
				tmp = posCen[k][i][dim];    // calculate the bin by dim;
				while (tmp > Lbox[dim]) 
				{ 
					tmp = tmp - Lbox[dim];
				}
				while (tmp < 0) 
				{ 
					tmp = tmp + Lbox[dim]; 
				}
				me = (tmp - lowPosition) / dbin;
				//me=(posCen[k][i][dim]-lowPosition)/dbin;
				if (me<nbin && me>-1)
				{
					bin[me] += 1.0;
					for (j = 0; j < halfIndex; j++)
					{
						r[0] = posCen[k + j][i][0] - posCen[k][i][0];
						r2 = r[0] * r[0];
						r[1] = posCen[k + j][i][1] - posCen[k][i][1];
						r2 = r[1] * r[1] + r2;
						msd[j][me] += r2;
					}
				}
			}
		}
	}
	if (region == 5)
	{
		//calculate the XYZ msd
		Dim_factor = 6;
		for (i = 0; i < nm; i++)
		{
			for (k = 0; k < halfIndex; k++)
			{
				tmp = posCen[k][i][dim];    // calculate the bin by dim;
				while (tmp > Lbox[dim]) 
				{ 
					tmp = tmp - Lbox[dim]; 
				}
				while (tmp < 0) 
				{ 
					tmp = tmp + Lbox[dim]; 
				}
				me = (tmp - lowPosition) / dbin;
				//me = (posCen[k][i][dim] - lowPosition)/dbin;
				if (me < nbin && me > -1)
				{
					bin[me] += 1.0;
					for (j = 0; j < halfIndex; j++)
					{
						r2 = 0;
						for (l = 0; l < 3; l++)
						{
							r[l] = periodicity(posCen[k + j][i][l] - posCen[k][i][l], Lbox[l]);
							r2 += r[l] * r[l];
						}
						msd[j][me] += r2;
					}
				}
			}
		}
	}
	for (i = 0; i < halfIndex; i++)
	{
		for (j = 0; j < nbin; j++)
		{
			if (bin[j] > 0.0)
			{
				msd[i][j] = msd[i][j] / bin[j];
			}
		}
	}
	return 0;
}

int load_positionCenter(t_trxframe fr, int nm, int atoms, int nStart, double** pos, double Lbox[3], double* deloc, double totNCM)
{
	int i, j, k, molN, m;
	double center[3], tmpPos, tmpCenter, halfLbox[3];
	for (k = 0; k < 3; k++)
		halfLbox[k] = Lbox[k] / 2;

	for (i = 0; i < nm; i++)
	{
		molN = nStart + i * atoms;

		for (k = 0; k < 3; k++) {
			tmpCenter = 0.0;  // initiate for each molecule
			for (j = 0; j < atoms; j++) {
				tmpPos = fr.x[molN + j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
				while (tmpPos - fr.x[molN][k] > halfLbox[k])
					tmpPos -= Lbox[k];
				while (tmpPos - fr.x[molN][k] < -halfLbox[k])
					tmpPos += Lbox[k];

				tmpCenter += tmpPos * deloc[j];
			}
			pos[i][k] = tmpCenter / totNCM;
		}
		// 	printf("i = %d pos = %f %f %f\n",i,pos[i][0],pos[i][1],pos[i][2]);
	}

	return 0;
}

int load_position(t_trxframe fr, int nm, int atoms, int nStart, double*** posmolecule, double Lbox[3])
{
	int i, j, k, molN;
	double center[3], tmpPos[3], halfLbox[3];
	for (k = 0; k < 3; k++)
		halfLbox[k] = Lbox[k] / 2;

	for (i = 0; i < nm; i++) {
		molN = nStart + i * atoms;
		for (j = 0; j < atoms; j++)
			for (k = 0; k < 3; k++) {
				tmpPos[k] = fr.x[molN + j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
				while (tmpPos[k] - fr.x[molN][k] > halfLbox[k])
					tmpPos[k] -= Lbox[k];
				while (tmpPos[k] - fr.x[molN][k] < -halfLbox[k])
					tmpPos[k] += Lbox[k];

				posmolecule[i][j][k] = tmpPos[k];
			}
	}

	return 0;
}

double periodicity(double dx, double box)
{
	while (dx > box / 2.0)
		dx -= box;
	while (dx < -box / 2.0)
		dx += box;
	return dx;
}

double* CreateVector(int rows)
{
	int	i;
	double* m;
	m = (double*) calloc((unsigned int)rows, sizeof(double*));
	return m;
}

double** CreateMatrix(int rows, int cols)
{
	int  i;
	double** m;
	m = (double**) calloc((unsigned int)rows, sizeof(double*));
	for (i = 0; i < rows; i++) {
		m[i] = (double*) calloc((unsigned int)cols, sizeof(double));
	}
	return m;
}
