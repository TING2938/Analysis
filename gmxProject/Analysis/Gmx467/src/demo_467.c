#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <gromacs/smalloc.h> 
#include <gromacs/copyrite.h>
#include <gromacs/statutil.h>
#include <gromacs/tpxio.h> 
#include <gromacs/macros.h>
#include <gromacs/index.h>    

int load_positionCenter(t_trxframe* fr, int nm, int atoms, int* index, double** posc, double* deloc);
int load_position(t_trxframe* fr, int nm, int atoms, int* index, double*** pos);
double periodicity(double dx, double box);
double* CreateVector(int cols);
double** CreateMatrix(int rows, int cols);
double*** CreateMatrix_3d(int m, int n, int p);
void getMolInfo(t_topology* top, int* ngx, int** index, int grp, int* nm, int* atoms, double** mass, double** charge);
void getLJParameter(t_topology* top, int** index, int grp1, int grp2, int atoms1, int atoms2, double*** c6, double*** c12);
void getMolInfo_diffType(t_topology* top, int* ngx, int** index, int grp, int** atomNum, int* nmols);
int load_positionCenter_diffType(t_trxframe* fr, t_topology* top, int nmols, int* eachAtom, int* index, double** posc, int ncm);

int main(int argc, char *argv[])
{
    const char *desc[] = {
      "this is a small test program meant to serve as a template "
    };

	int nbin = 100;
	real lowPos = 0.0, upPos = 30.0;

	t_pargs pa[] = {
	  { "-nbin",FALSE, etINT, {&nbin},
		"number of bin sets for number"
	  },
	  { "-low",FALSE, etREAL, {&lowPos},
		"low position of region of molecule/ion (nm)"
	  },
	  { "-up",FALSE, etREAL, {&upPos},
		"up position of region of molecule/ion (nm)"
	  }
	};

    t_filenm           fnm[] = {
        { efTRX, "-f", NULL,  ffREAD },
        { efTPX, NULL, NULL,  ffREAD },
        { efNDX, NULL, NULL,  ffOPTRD },
        { efXVG, NULL, "density", ffWRITE }
    };

    /* GMX-467 support data structure */
    t_topology*     top;
    int             ePBC;   
    output_env_t    oenv;
    t_trxframe      fr;
    t_trxstatus*    status;
    int             flags = TRX_READ_X;   /* read only position     */
    char**          grpname;              /* group names            */
    int             ngrps = 2;            /* nr. of group(s)        */
    int**           index;                /* indices for all groups */
    int*            ngx;                  /* sizes of groups        */

#define NFILE asize(fnm)
    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    top = read_top(ftp2fn(efTPX, NFILE, fnm), &ePBC); /* read topology file */
    snew(grpname, ngrps);
    snew(index, ngrps);
    snew(ngx, ngrps);
	printf("\nSpecify %d group%s to analysis:\n", ngrps, (ngrps > 1) ? "s" : "");
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname); 
	/* ======================= Main body of code ========================= */

	int nm1, nm2, atoms1, atoms2;
	// double ***pos1, ***pos2;                         /* size: [nm, atoms, 3]     */
	// double **posc1, **posc2;                         /* size: [nm, 3]            */
	// double *mass1, *mass2, *charge1, *charge2;       /* size: [atoms]            */
	// double **c6, **c12;                           /* LJ parameter. size: [atoms1, atoms2] */

	// getMolInfo(top, ngx, index, 0, &nm1, &atoms1, &mass1, &charge1);
	// getMolInfo(top, ngx, index, 1, &nm2, &atoms2, &mass2, &charge2);
	// getLJParameter(top, index, 0, 1, atoms1, atoms2, &c6, &c12);

    /* if the selection contains molecules with different type, please use function below to get 
     * number of molecules (`nmols`) and atoms of each molecule (`eachAtom`), in this case, you 
     * should comment the above codes that contain `getMolInfo` and `getLJParameter`;
    */
    
    int* eachAtom;      // nr. of atoms of each molecule
    int nmols;          // nr. of molecules in the selection
    getMolInfo_diffType(top, ngx, index, 0, &eachAtom, &nmols);
    double** posc = CreateMatrix(nmols, 3);

    // load_positionCenter_diffType(&fr, top, nmols, eachAtom, index[0], posc, 2);
    
    

	read_first_frame(oenv, &status, ftp2fn_null(efTRX, NFILE, fnm), &fr, flags);

    // pos1 = CreateMatrix_3d(nm1, atoms1, 3);
    // pos2 = CreateMatrix_3d(nm2, atoms2, 3); 
    // posc1 = CreateMatrix(nm1, 3);   
    // posc2 = CreateMatrix(nm2, 3);

    do {
        load_positionCenter_diffType(&fr, top, nmols, eachAtom, index[0], posc, 2);
		// load_position(&fr, nm1, atoms1, index[0], pos1);
		// load_position(&fr, nm2, atoms2, index[1], pos2);
        // load_positionCenter(&fr, nm1, atoms1, index[0], posc1, mass1);
		// load_positionCenter(&fr, nm2, atoms2, index[1], posc2, mass2);
		
    } while (read_next_frame(oenv, status, &fr));

    
    thanx(stderr);
    return 0;
}


//reference function
int load_positionCenter(t_trxframe* fr, int nm, int atoms, int *index, double **pos, double *deloc)
{
    int i, j, k, molN;
    double tmpPos, tmpCenter, halfLbox[3];
	double totNCM = 0.0;
	for (i = 0; i != atoms; ++i) {
		totNCM += deloc[i];
	}

    for (k = 0; k < 3; k++)
        halfLbox[k] = fr->box[k][k] / 2;

    for (i = 0; i < nm; i++) {
        molN = index[i * atoms];

        for (k = 0; k < 3; k++) {
            tmpCenter = 0.0;  // initiate for each molecule
            for (j = 0; j < atoms; j++) {
                tmpPos = fr->x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos - fr->x[molN][k] > halfLbox[k])
                    tmpPos -= fr->box[k][k];
                while (tmpPos - fr->x[molN][k] < -halfLbox[k])
                    tmpPos += fr->box[k][k];

                tmpCenter += tmpPos * deloc[j];
            }
            pos[i][k] = tmpCenter / totNCM;
        }
    }

    return 0;
}

int load_position(t_trxframe* fr, int nm, int atoms, int *index, double ***posmolecule)
{
    int i, j, k, molN;
    double tmpPos[3], halfLbox[3];
    for (k = 0; k < 3; k++)
        halfLbox[k] = fr->box[k][k] / 2;

    for (i = 0; i < nm; i++) {
        molN = index[i * atoms];
        for (j = 0; j < atoms; j++)
            for (k = 0; k < 3; k++) {
                tmpPos[k] = fr->x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos[k] - fr->x[molN][k] > halfLbox[k])
                    tmpPos[k] -= fr->box[k][k];
                while (tmpPos[k] - fr->x[molN][k] < -halfLbox[k])
                    tmpPos[k] += fr->box[k][k];

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
    return (double*)malloc(rows * sizeof(double));
}

double** CreateMatrix(int rows, int cols)
{
    double  **m  = (double**)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        m[i] = (double*)malloc(cols * sizeof(double));
    }
    return m;
}

double*** CreateMatrix_3d(int m, int n, int p)
{
	double*** m2 = (double***)malloc(m * sizeof(double**));
	for (int i = 0; i < m; i++)
		m2[i] = CreateMatrix(n, p);
	return m2;
}

void getMolInfo(t_topology* top, int* ngx, int** index, int grp, int* nm, int* atoms, double** mass, double** charge)
{
	int i;
    *atoms = 0;
    for (;;)
    {
        if (top->atoms.atom[index[grp][*atoms]].resind == top->atoms.atom[index[grp][0]].resind)
            (*atoms)++;
        else
            break;
    }
	*nm = ngx[grp] / *atoms;

	printf("\n### Group %d: \n", grp);
	printf("# Nr. of atoms in molecule : %d\n", *atoms);
	printf("# Nr. of molecules in group: %d\n", *nm); 
	printf("# nr.   name     mass   charge\n");
	printf("------------------------------\n");
	double sumMass = 0.0, sumCharge = 0.0;
	snew(*mass, *atoms);
	snew(*charge, *atoms);
	for (i = 0; i != *atoms; ++i) {
		(*mass)[i] = top->atoms.atom[index[grp][i]].m;
		(*charge)[i] = top->atoms.atom[index[grp][i]].q;
		sumMass += (*mass)[i];
		sumCharge += (*charge)[i];
		printf("%5d%7s%9.3f%9.3f\n",i, *top->atoms.atomname[index[grp][i]], (*mass)[i], (*charge)[i]);
	}
	printf("------------------------------\n");
	printf("%12s%9.3f%9.3f\n", "Total", sumMass, sumCharge);

}

void getLJParameter(t_topology* top, int** index, int grp1, int grp2, int atoms1, int atoms2, double*** c6, double*** c12)
{
	int n1, n2, i, j;
	int ntypes = top->atomtypes.nr;
	int type;
	*c6 = CreateMatrix(atoms1, atoms2);
	*c12 = CreateMatrix(atoms1, atoms2); 
	printf("\n### LJ Parameters(Group %d and %d):\n", grp1, grp2);
	printf("[  at1  at2]: %12s %12s\n", "c6", "c12");
	printf("---------------------------------------\n");
	int n = 0;
	for (i = 0; i != atoms1; ++i) {
		for (j = 0; j != atoms2; ++j) {
			n1 = index[grp1][i];
			n2 = index[grp2][j];
			type = ntypes * (top->atoms.atom[n1].type) + (top->atoms.atom[n2].type);
			(*c6)[i][j] = top->idef.iparams[type].lj.c6;
			(*c12)[i][j] = top->idef.iparams[type].lj.c12;
			if (n++ <= 12) {
				printf("[%5s%5s]: %12.4e %12.4e\n", *top->atoms.atomname[n1], *top->atoms.atomname[n2], (*c6)[i][j], (*c12)[i][j]);
			}
		}
	}
	if (atoms1 * atoms2 > 12) {
		printf("(.........)");
	}
	printf("\n\n");
}

void getMolInfo_diffType(t_topology* top, int* ngx, int** index, int grp, int** atomNum, int* nmols)
{
    int* tmpAtoms = (int*)malloc(sizeof(int) * ngx[grp]);
    for (int i = 0; i < ngx[grp]; i++)
    {
        tmpAtoms[i] = 0;
    }

    int firstResi;
    int j = 0;
    firstResi = top->atoms.atom[index[grp][0]].resind;

    for (int i = 0; i < ngx[grp]; i++)
    {
        if (top->atoms.atom[index[grp][i]].resind == firstResi)
        {
            tmpAtoms[j]++;
        }
        else
        {
            firstResi = top->atoms.atom[index[grp][i]].resind;
            j++;
            tmpAtoms[j] = 1;
        } 
    }
    *nmols = j + 1;
    *atomNum = (int*)malloc(sizeof(int) * (*nmols));
    for (int i = 0; i != *nmols; i++)
    {
        (*atomNum)[i] = tmpAtoms[i];
    }
    free(tmpAtoms);
    printf("\nThere are %d molecules in your selection %d\n", *nmols, grp);
}

int load_positionCenter_diffType(t_trxframe* fr, t_topology* top, int nmols, int* eachAtom, int* index, double** posc, int ncm)
{
    static gmx_bool bFirst = 1;
    static double** deloc;
    static double* totNCM;
    if (bFirst)
    {
        bFirst = 0; // just do it once
        totNCM = (double*)malloc(sizeof(double) * nmols);
        deloc = (double**)malloc(sizeof(double*) * nmols);
        int atomIndex = 0;
        for (int i = 0; i < nmols; i++)
        {
            deloc[i] = (double*)malloc(sizeof(double) * eachAtom[i]);
            totNCM[i] = 0;

            if (ncm == 0) // number center
            {
                for (int j = 0; j < eachAtom[i]; j++)
                {
                    deloc[i][j] = 1;
                }
                totNCM[i] = eachAtom[i];
            }
            else if (ncm == 1) // charge center
            {
                for (int j = 0; j < eachAtom[i]; j++)
                {
                    deloc[i][j] = top->atoms.atom[index[atomIndex++]].q;
                    totNCM[i] += deloc[i][j];
                }
            }
            else if (ncm == 2) // mass center
            {
                for (int j = 0; j < eachAtom[i]; j++)
                {
                    deloc[i][j] = top->atoms.atom[index[atomIndex++]].m;
                    totNCM[i] += deloc[i][j];
                }
            }
            else
            {
                printf("Error of selected ncm, pleace check your commandline argument!\n");
                exit(-1);
            }
        }
    }

    int i, j, k, molN;
    double tmpPos, tmpCenter, halfLbox[3];
    int firstIndex = 0;

    for (k = 0; k < 3; k++)
        halfLbox[k] = fr->box[k][k] / 2;

    for (i = 0; i < nmols; i++)
    {
        molN = index[firstIndex];

        for (k = 0; k < 3; k++)
        {
            tmpCenter = 0.0;  // initiate for each molecule
            for (j = 0; j < eachAtom[i]; j++)
            {
                tmpPos = fr->x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos - fr->x[molN][k] > halfLbox[k])
                    tmpPos -= fr->box[k][k];
                while (tmpPos - fr->x[molN][k] < -halfLbox[k])
                    tmpPos += fr->box[k][k];
                tmpCenter += tmpPos * deloc[i][j];
            }
            posc[i][k] = tmpCenter / totNCM[i];
        }
        firstIndex += eachAtom[i];
    }
    return 0;
}

