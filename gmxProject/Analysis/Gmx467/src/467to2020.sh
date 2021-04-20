#!/usr/bin/env bash
### convert GMX467 trajectory analysis code to GMX2020
## yeting
## 2020/09/01
## yeting2938@hust.edu.cn

# backup
cp $1 ${1}".bak"

fnm=$1
# for new header
sed -i '/#include/d' $fnm
sed -i '1i\#include <gromacs/commandline/pargs.h>' $fnm
sed -i '1i\#include <gromacs/fileio/trxio.h>' $fnm
sed -i '1i\#include <gromacs/fileio/tpxio.h>' $fnm
sed -i '1i\#include <gromacs/mdtypes/inputrec.h>' $fnm
sed -i '1i\#include <gromacs/utility/fatalerror.h>' $fnm
sed -i '1i\#include <gromacs/topology/index.h>' $fnm
sed -i '1i\#include <gromacs/topology/topology.h>' $fnm
sed -i '1i\#include <gromacs/utility/smalloc.h>' $fnm
sed -i '1i\#include <gromacs/trajectory/trajectoryframe.h>' $fnm
sed -i '1i\#include <gromacs/commandline/cmdlineinit.h>' $fnm
sed -i '1i\#include <gromacs/utility/arraysize.h>' $fnm
sed -i '1i\#include <gromacs/fileio/xvgr.h>' $fnm
sed -i '1i\#include <gromacs/fileio/confio.h>' $fnm

# change
sed -i 's/output_env_t/gmx_output_env_t*/g' $fnm
sed -i 's/atom_id/int/g' $fnm
sed -i '/CopyRight(/d' $fnm
sed -i '/thanx(/d' $fnm
sed -i 's/main(/main_func(/g' $fnm
sed -i 's/PCA_BE_NICE/PCA_CAN_TIME/g' $fnm
sed -i 's/efTPX/efTPR/g' $fnm

line=`grep -A 10 'parse_common_args' $fnm | /bin/grep ';' | head -n 1 | xargs -n 1 -i /bin/grep -n {} $fnm`
sed -i "${line%%:*}s/;/)\n{\n\texit(0);\n}/g" $fnm
sed -i 's/parse_common_args/if(!parse_common_args/g' $fnm
echo -e '\nint main(int argc, char** argv)\n{\n\treturn gmx_run_cmain(argc, argv, &main_func);\n}\n' >> $fnm

