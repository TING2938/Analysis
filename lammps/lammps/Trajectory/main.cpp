#include "trajectory.h"

int main(int argc, char** argv)
{
	Getopt opt(argc, argv);
	int aa = 0;
	opt(aa, "-a");

	Trajectory trj(argc, argv);
	trj.readFirstFrame();
	do {

	} while (trj.readNextFrame());
	

}

