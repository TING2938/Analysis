#include "src/conp.h"

int main(int argc, char** argv)
{
	itp::Getopt subprogram(argc, argv);

	Conp conp(argc, argv);

	subprogram.addSubProgram("matrix", "Calculate Matrix for ECPM", 
		[&conp] {conp.calc_Matrix(); }
		);

	subprogram.addSubProgram("potFile", "Get \"getPot_parameters.dat\" for different voltage", 
		[&conp] { conp.get_getPotFile(); }
		);

	subprogram.finish();
}

