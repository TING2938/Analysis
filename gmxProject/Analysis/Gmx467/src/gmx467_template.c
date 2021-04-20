#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <ctype.h>
#include <gromacs/sysstuff.h>
#include <string.h>
#include <gromacs/string2.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/macros.h>
#include <gromacs/vec.h>
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


int main(int argc, char *argv[])
{
    const char        *desc[] = {
        "do some analysis about trajectory.",
    };

    output_env_t       oenv;
    char           **grpname;      /* groupnames               */
    int                 *ngx;      /* sizes of groups          */
    t_topology          *top;      /* topology                 */
    int                 ePBC;
    atom_id          **index;      /* indices for all groups   */
    t_trxframe            fr;
    t_trxstatus      *status;
    FILE *fp = NULL;

    int   flags = TRX_READ_X;  /* flag . which (x, v, or f) to read*/
    static int         ngrps = 2;  /* nr. of groups */
    static gmx_bool    bCenter = FALSE;

    t_pargs            pa[] = {
        { "-center",  FALSE, etBOOL, {&bCenter},
          "Shift the center of mass along the axis to zero."}
    };
    t_filenm           fnm[] = { /* files io  */
        { efTRX, "-f", NULL,  ffREAD },
        { efNDX, NULL, NULL,  ffOPTRD },
        { efTPX, NULL, NULL,  ffREAD },
        { efXVG, "-o", "density", ffWRITE },
    };
#define NFILE asize(fnm)
    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    top = read_top(ftp2fn(efTPX, NFILE, fnm), &ePBC); /* read topology file */
    snew(grpname, ngrps);
    snew(index, ngrps);
    snew(ngx, ngrps);
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);
    read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);
    if (opt2bSet("-o", NFILE, fnm)) {
        fp = xvgropen(opt2fn("-o", NFILE, fnm), "Average distance",
                      "Time [ps]", "Distance [nm]", oenv);
    }
    /****************** do some analysis     *****************************************/
    do {
        double xx = 0;
        for (size_t n = 0; n != ngrps; n++) /* loop over all groups */
        {
            for (size_t i = 0; i != ngx[n]; i++) /* loop over all atoms in index file */
            {
                atom_id ndx = index[n][i];  /* index of each atom */
                //xx = fr.x[ndx][0];
                //printf("%5d", ndx);
                // auto vv = fr.v[ndx];
                // auto ff = fr.f[ndx];
            }
        }
        fprintf(fp, "%10.3f", fr.time);
        fprintf(fp, "%10.3f\n", fr.x[1][1]);
    } while (read_next_frame(oenv, status, &fr));
    /*********************************************************************************/
    if (fp) { ffclose(fp); }
    thanx(stderr);
    return 0;
}
