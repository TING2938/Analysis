#ifdef HAVE_CONFIG_H
#include <gromacs/config.h>
#endif
#include <math.h>
#include <ctype.h>

#include <gromacs/sysstuff.h>
#include <string.h>
#include <gromacs/string2.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/macros.h>
#include <gromacs/gstat.h>
#include <gromacs/vec.h>
#include <gromacs/xvgr.h>
#include <gromacs/pbc.h>
#include <gromacs/copyrite.h>
#include <gromacs/futil.h>
#include <gromacs/statutil.h>
#include <gromacs/index.h>
#include <gromacs/tpxio.h>
#include <gromacs/physics.h>
#include <gromacs/gmx_ana.h>

/****************************************************************************/
/* This program calculates the diffusion by msd.                  */
/* fengjie peng july 2016                                                 */
/****************************************************************************/

static int  dim = 2;
static real zmin = 0;
static real zmax = 0;
static real c1 = 0;
static real c2 = 0;
static real Cr = 0;
static real CR = 0.6;
static real rcut = 1.0;
static int  maxfr = 100000;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static int  rmind = 0;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static real kB = 1.3806505;/* Boltzmann constant without 10^-23 */
static real NA = 6.02214129;  /* Avogadro constant without 10^23 */
static real fcomb = 138.935485; /* kJ mol^−1 nm e^−2 factor in coulomb interaction */
static const char* normtype[] = { NULL,"no","x","y","z",NULL };
static const char* axtitle[] = { NULL,"no","x","y","z",NULL };

double* CreateVector(int rows);
double** CreateMatrix(int rows, int cols);
double*** CreateMatrix_3d(int m, int n, int p);
double pbc(double dx, double box);
int msd_compute(int nm, int step, int nbin, double** bin, double*** posCen, double** msd, double Lbox[3], double dbin, double lowPosition);
int load_positionCenter(rvec* x0, int nm, int napm, int nStart, double** pos, double Lbox[3], double* deloc, double totMNC);
int load_position(rvec* x0, int nm, int napm, int nStart, double*** posmolecule, double Lbox[3]);
static void corr_print(const char* fn, const char* title, const char* yaxis, double dt, int step, int nbin, double** msd,
	char* grpname[], const output_env_t oenv);
void do_corr(const char* trx_file, const char* ndx_file, const char* msd_file,
	int* tot_atom, char* grpname[], gmx_bool bMOL, t_topology* top, int nbin, const output_env_t oenv)
{
	rvec* x0;      /* coordinates without pbc */
	matrix       box;     /* box (3x3) */
	t_trxstatus* status;
	real        t, tpf;
	int         natoms;  /* nr. atoms in trj */
	int         i, j, k, k1, k2, step, nm, nstart, napm;
	// int         flags = TRX_READ_X;   /* read only position */
	double      Lbox[3], dt, halfLbox[3];
	atom_id** index;
	double** msd, *** posCen, ** rmPos, ** tempPos, dbin, * mass, tmass, ** numbin;                /* some 2d-matrix to be used */

	// read_first_frame(oenv,&status,trx_file,&fr,flags);
	natoms = read_first_x(oenv, &status, trx_file, &t, &x0, box);
	Lbox[0] = box[XX][XX]; halfLbox[0] = 0.5 * Lbox[0];
	Lbox[1] = box[YY][YY]; halfLbox[1] = 0.5 * Lbox[1];
	Lbox[2] = box[ZZ][ZZ]; halfLbox[2] = 0.5 * Lbox[2];
	tpf = t;

	if (zmax == 0) { zmax = Lbox[dim]; }
	fprintf(stdout, " \n ");
	snew(index, 1);
	get_index(&top->atoms, ndx_file, 1, tot_atom, index, grpname);
	nstart = index[0][0];
	napm = top->mols.index[1 + top->atoms.atom[nstart].resind] - top->mols.index[top->atoms.atom[nstart].resind];
	if (!bMOL) { napm = 1; }
	nm = tot_atom[0] / napm;
	dbin = (zmax - zmin) / nbin;
	msd = CreateMatrix(maxfr, nbin);
	numbin = CreateMatrix(maxfr, nbin);
	rmPos = CreateMatrix(nm, 3);
	tempPos = CreateMatrix(nm, 3);
	posCen = CreateMatrix_3d(maxfr, nm, 3);
	mass = CreateVector(napm);
	tmass = 0.0;
	for (i = 0; i < napm; i++) {
		mass[i] = top->atoms.atom[nstart + i].m;
		tmass += mass[i];
	}
	for (i = 0; i < nbin; i++) {
		for (j = 0; j < maxfr; j++) {
			msd[j][i] = 0.0;
			numbin[j][i] = 0.0;
		}
	}
	for (i = 0; i < nm; i++) {
		for (j = 0; j < 3; j++) {
			rmPos[i][j] = 0.0;
		}
	}

	step = 0;
	do {

		load_positionCenter(x0, nm, napm, nstart - rmind, posCen[step], Lbox, mass, tmass);
		for (i = 0; i < nm; i++) {
			for (j = 0; j < 3; j++) {
				if (step == 0) {
					tempPos[i][j] = posCen[step][i][j];
				}
				if (step > 0) {
					if ((posCen[step][i][j] - tempPos[i][j]) < -halfLbox[j])
						rmPos[i][j] += Lbox[j];
					if ((posCen[step][i][j] - tempPos[i][j]) > halfLbox[j])
						rmPos[i][j] -= Lbox[j];
					tempPos[i][j] = posCen[step][i][j];
					posCen[step][i][j] += rmPos[i][j];
				}
			}
		}
		// msd_compute(nm,step,nbin,numbin,posCen,msd,Lbox,dbin,zmin);
		step++;
		dt = t - tpf;
		tpf = t;
	} while (read_next_x(oenv, status, &t, natoms, x0, box));
	msd_compute(nm, step, nbin, numbin, posCen, msd, Lbox, dbin, zmin);
	for (i = 0; i < nbin; i++) {
		for (j = 0; j < step / 2; j++) {
			if (numbin[j][i] > 10.0) {
				msd[j][i] /= numbin[j][i];
			}
			else { msd[j][i] = 0; }
		}
	}

	// fprintf(stdout," \n");
	corr_print(msd_file, "Mean Square Displacement.", "MSD (nm\\S2\\N).", dt, step, nbin, msd, grpname, oenv);
}

int main(int argc, char* argv[])
{
	const char* desc[] = {
	  "This program calculates the diffusion by msd.",
	};

	static int  nbin = 10;
	// static real beginfit   = -1; 
	// static real endfit     = -1; 
	static gmx_bool bMOL = TRUE;
	// static gmx_bool bRmCOMM    = FALSE;
	t_pargs pa[] = {
	  { "-mol",    FALSE, etBOOL, {&bMOL},
		"count by molecules or atoms, default is treated as molecules" },
	  { "-type",    FALSE, etENUM, {normtype},
		"Compute diffusion coefficient in one direction" },
	  { "-lateral", FALSE, etENUM, {axtitle},
		"Calculate the lateral diffusion in a plane perpendicular to" },
	  { "-rm",FALSE, etINT, {&rmind},
		"number of atoms didn't storage in xtc file before selected group" },
	  { "-maxfr",FALSE, etINT, {&maxfr},
		"maximum of frames in trajectory" },
	  { "-sl",  FALSE, etINT, {&nbin},
		"Divide the box in #nr slices." },
	  { "-minz",FALSE, etREAL, {&zmin},
		"min z to calculate (nm)"},
	  { "-maxz",FALSE, etREAL, {&zmax},
		"max z to calculate (nm)"},
		{ "-c1",FALSE, etREAL, {&c1},
		"center 1 to calculate (nm)"},
		{ "-c2",FALSE, etREAL, {&c2},
		"center 2 to calculate (nm)"},
		{ "-Cr",FALSE, etREAL, {&Cr},
		"R min to calculate (nm)"},
		{ "-CR",FALSE, etREAL, {&CR},
		"R max to calculate (nm)"},
	  { "-dim",FALSE, etINT, {&dim},
		"dimension to calculate, x(0), y(1), z(2)" },
		// { "-beginfit",FALSE, etTIME, {&beginfit},
		  // "Start time for fitting the MSD (%t), -1 is 10%" },
		// { "-endfit",FALSE, etTIME, {&endfit},
		  // "End time for fitting the MSD (%t), -1 is 90%" }
	};

	t_filenm fnm[] = {
	  { efTRX, NULL, "run",  ffREAD },
	  { efTPS, NULL, "run",  ffREAD },
	  { efNDX, NULL, NULL,  ffOPTRD },
	  { efXVG, NULL, "msd", ffWRITE },
	};
#define NFILE asize(fnm)

	t_topology* top;
	int         ePBC;
	const char* trx_file, * tps_file, * ndx_file, * msd_file;
	int* gnx; /* the selected groups' sizes */
	char** grpname;
	output_env_t oenv;

	CopyRight(stderr, argv[0]);

	parse_common_args(&argc, argv,
		PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_BE_NICE,
		NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);
	trx_file = ftp2fn_null(efTRX, NFILE, fnm);
	tps_file = ftp2fn_null(efTPS, NFILE, fnm);
	ndx_file = ftp2fn_null(efNDX, NFILE, fnm);
	msd_file = ftp2fn_null(efXVG, NFILE, fnm);


	snew(gnx, 1);
	snew(grpname, 1);
	top = read_top(tps_file, &ePBC);     /* read topology file */
   // get_index(&top->atoms,ndx_file,1,gnx,index,grpname);
  /* the function you are going to use */
	do_corr(trx_file, ndx_file, msd_file, gnx, grpname, bMOL,
		top, nbin, oenv);
	/*  the function end  */
	thanx(stderr);

	return 0;
}

double* CreateVector(int rows)
{
	double* m;

	m = calloc((unsigned int)rows, sizeof(double));

	return m;
}
double** CreateMatrix(int rows, int cols)
{
	int  i;
	double** m;

	m = calloc((unsigned int)rows, sizeof(double*));
	for (i = 0; i < rows; i++) {
		m[i] = calloc((unsigned int)cols, sizeof(double));
	}

	return m;
}
double*** CreateMatrix_3d(int m, int n, int p)
{
	int i;
	double*** m2;
	m2 = calloc((unsigned int)m, sizeof(double**));
	for (i = 0; i < m; i++)
		m2[i] = CreateMatrix(n, p);
	return m2;
}
double pbc(double dx, double box)
{
	while (dx > box / 2.0)
		dx -= box;
	while (dx < -box / 2.0)
		dx += box;

	return dx;
}
int load_positionCenter(rvec* x0, int nm, int napm, int nStart, double** pos, double Lbox[3], double* deloc, double totMNC)
{
	int i, j, k, molN, m;
	double center[3], tmpPos, tmpCenter, halfLbox[3];
	for (k = 0; k < 3; k++)
		halfLbox[k] = Lbox[k] / 2;

	for (i = 0; i < nm; i++) {
		molN = nStart + i * napm;

		for (k = 0; k < 3; k++) {
			tmpCenter = 0.0;  // initiate for each molecule
			for (j = 0; j < napm; j++) {
				tmpPos = x0[molN + j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
				while (tmpPos - x0[molN][k] > halfLbox[k])
					tmpPos -= Lbox[k];
				while (tmpPos - x0[molN][k] < -halfLbox[k])
					tmpPos += Lbox[k];

				tmpCenter += tmpPos * deloc[j];
			}
			pos[i][k] = tmpCenter / totMNC;
			// move molecules into box.
			while (pos[i][k] > Lbox[k])
				pos[i][k] -= Lbox[k];
			while (pos[i][k] < 0)
				pos[i][k] += Lbox[k];
		}
	}

	return 0;
}
int load_position(rvec* x0, int nm, int napm, int nStart, double*** posmolecule, double Lbox[3])
{
	int i, j, k, molN;
	double center[3], tmpPos[3], halfLbox[3];
	for (k = 0; k < 3; k++)
		halfLbox[k] = Lbox[k] / 2;

	for (i = 0; i < nm; i++) {
		molN = nStart + i * napm;
		for (j = 0; j < napm; j++)
			for (k = 0; k < 3; k++) {
				tmpPos[k] = x0[molN + j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
				while (tmpPos[k] - x0[molN][k] > halfLbox[k])
					tmpPos[k] -= Lbox[k];
				while (tmpPos[k] - x0[molN][k] < -halfLbox[k])
					tmpPos[k] += Lbox[k];

				posmolecule[i][j][k] = tmpPos[k];
			}
	}

	return 0;
}


int msd_compute(int nm, int step, int nbin, double** bin, double*** posCen, double** msd, double Lbox[3], double dbin, double lowPosition)
{
	int i, k1, k2, m, me, me1, totbin, conut;
	double r, r2;
	double tmpPos1, tmpPos2, Rsquare;

	totbin = Lbox[dim] / dbin;
	for (i = 0; i < nm; i++) {
		for (k1 = 0; k1 < step / 2; k1++) {
			me = (int)((posCen[k1][i][dim] - lowPosition) / dbin + 1) - 1;
			// me =  floor((posCen[k1][i][2]-lowPosition)/dbin+0.5);
			// by +1 and -1 to avoid (-1,0) -> 0 , another way is using floor() function;
			conut = 0;
			// if(step ==1)
			// fprintf(stdout,"  me  posCen[step-k][i][2] dbin %4i %4f %4f \n",me,posCen[step-k][i][2],dbin);
			while (me < 0)
				me += totbin;
			while (me >= totbin)
				me -= totbin;

			tmpPos1 = posCen[k1][i][XX] - c1;
			tmpPos2 = posCen[k1][i][YY] - c2;
			Rsquare = tmpPos1 * tmpPos1 + tmpPos2 * tmpPos2;
			if (me<nbin && me>-1 && Cr*Cr <= Rsquare && Rsquare <= CR*CR) {
				for (k2 = 0; k2 < step / 2; k2++) {
					me1 = (int)((posCen[k1 + k2][i][dim] - lowPosition) / dbin + 1) - 1;
					// me1 = floor((posCen[k1+k2][i][2]-lowPosition)/dbin+0.5);
					while (me1 < 0)
						me1 += totbin;
					while (me1 >= totbin)
						me1 -= totbin;
					if (me1 == me) { conut++; }
					if (conut == k2 + 1) {
						bin[k2][me] += 1.0;
						if (normtype[0][0] == 'n' && axtitle[0][0] == 'n') {
							for (m = 0; m < 3; m++) {
								r = posCen[k1 + k2][i][m] - posCen[k1][i][m];
								r2 = r * r;
								// fprintf(stdout,"  r2 %4f ",r2);
								msd[k2][me] += r2;
							}
						}

						if (normtype[0][0] == 'x' && axtitle[0][0] == 'n') {
							r = posCen[k1 + k2][i][0] - posCen[k1][i][0];
							r2 = r * r;
							msd[k2][me] += r2;
						}

						if (normtype[0][0] == 'y' && axtitle[0][0] == 'n') {
							r = posCen[k1 + k2][i][1] - posCen[k1][i][1];
							r2 = r * r;
							msd[k2][me] += r2;
						}

						if (normtype[0][0] == 'z' && axtitle[0][0] == 'n') {
							r = posCen[k1 + k2][i][2] - posCen[k1][i][2];
							r2 = r * r;
							msd[k2][me] += r2;
						}

						if (normtype[0][0] == 'n' && axtitle[0][0] == 'x') {
							r = posCen[k1 + k2][i][1] - posCen[k1][i][1];
							r2 = r * r;
							msd[k2][me] += r2;
							r = posCen[k1 + k2][i][2] - posCen[k1][i][2];
							r2 = r * r;
							msd[k2][me] += r2;
						}

						if (normtype[0][0] == 'n' && axtitle[0][0] == 'y') {
							r = posCen[k1 + k2][i][0] - posCen[k1][i][0];
							r2 = r * r;
							msd[k2][me] += r2;
							r = posCen[k1 + k2][i][2] - posCen[k1][i][2];
							r2 = r * r;
							msd[k2][me] += r2;
						}

						if (normtype[0][0] == 'n' && axtitle[0][0] == 'z') {
							r = posCen[k1 + k2][i][1] - posCen[k1][i][1];
							r2 = r * r;
							msd[k2][me] += r2;
							r = posCen[k1 + k2][i][0] - posCen[k1][i][0];
							r2 = r * r;
							msd[k2][me] += r2;
						}
					}
					else k2 = step / 2;
				}
			}
		}
	}
	return 0;
}

static void corr_print(const char* fn, const char* title, const char* yaxis, double dt, int step, int nbin, double** msd,
	char* grpname[], const output_env_t oenv)
{
	FILE* out;
	int  i, j;

	out = xvgropen(fn, title, output_env_get_xvgr_tlabel(oenv), yaxis, oenv);
	// fprintf(out,"# MSD gathered over %s with 100 restarts\n",
		// output_env_get_time_unit(oenv));
	// fprintf(out,"# Diffusion constants fitted from time %g to %g %s\n",
		// beginfit,endfit,output_env_get_time_unit(oenv));
	for (i = 0; i < step / 2; i++) {
		fprintf(out, "  %4.1f    ", i * dt);
		for (j = 0; j < nbin; j++) {
			fprintf(out, "%10g   ", msd[i][j]);
		}
		fprintf(out, "\n");
	}
	ffclose(out);
}


