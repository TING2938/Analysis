#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "core_getopt.h"

class Trajectory
{
public:
	Trajectory(int argc, char** argv);
	~Trajectory();

	void readFirstFrame();
	bool readNextFrame();

	double* getData(std::string sign);
	int getTime();

private:
	double                      **data;
	int                         natoms;
	int                         ncolumn; 
	int                         curTime;
	double                      Lbox[3][2];
	std::vector<std::string>    headline;
	std::ifstream               fid; 
};

