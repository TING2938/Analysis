#include "trajectory.h" 

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

Trajectory::Trajectory(int argc, char** argv)
{ 
	Getopt opt(argc, argv);
	std::string fileName;
	opt(fileName, "-f");
	fid.open(fileName);
	if (!fid.is_open()) {
		std::cout << "Cannot open file(\"" << fileName << "\")!\nExit" << std::endl;
		exit(-1);
	}

	printf("# ----------------------- Lammps Trajectory Analysis ----------------------- #\n");
	
	char pwd[FILENAME_MAX];
	GetCurrentDir(pwd, sizeof(pwd));
	printf("### Working dir : %s\n", pwd);
	printf("### Command line: %s\n", opt.toString().c_str());
	printf("\n");
}

void Trajectory::readFirstFrame()
{
	std::string tmp, line;
	std::stringstream ss;

	std::getline(fid, line);

	// get current time;
	std::getline(fid, line);
	ss.str(line);
	ss >> curTime;
	ss.clear();

	std::getline(fid, line);

	// get nr. atoms;
	std::getline(fid, line);
	ss.str(line);
	ss >> natoms;
	ss.clear();

	if (natoms == 0) {
		std::cout << "No atoms!\nExit" << std::endl;
		exit(-1);
	}

	std::getline(fid, line);

	// get box size;
	for (int i = 0; i != 3; ++i) {
		std::getline(fid, line);
		ss.str(line);
		ss >> Lbox[i][0] >> Lbox[i][1];
		ss.clear();
	}

	// get headline;
	std::getline(fid, line);
	ss.str(line);
	ss >> tmp >> tmp;
	while (ss >> tmp)
	{
		headline.push_back(tmp);
	}
	ss.clear();
	ncolumn = (int) headline.size();

	if (ncolumn == 0) {
		std::cout << "No data!\nExit" << std::endl;
		exit(-1);
	}

	data = new double* [ncolumn];
	for (int i = 0; i != ncolumn; ++i) {
		data[i] = new double[natoms];
	}

	// get data;
	for (int i = 0; i != natoms; ++i) {
		std::getline(fid, line);
		ss.str(line);
		for (int j = 0; j != ncolumn; ++j) {
			ss >> data[j][i];
		}
		ss.clear();
	}

	printf("### Read First Frame\n");
	printf("Natoms     : %d\n", natoms);
	printf("Ncolumn    : %d\n", ncolumn);
	printf("Head Line  : ");
	for (auto&& i : headline)
		printf("%s ", i.c_str());
	printf("\n\n");
	printf("### Start Analysis\n");

}

bool Trajectory::readNextFrame()
{
	std::string tmp, line;
	std::stringstream ss;

	if (std::getline(fid, line)) {

		// get current time;
		std::getline(fid, line);
		ss.str(line);
		ss >> curTime;
		ss.clear();

		for (int i = 0; i != 3; ++i)
			std::getline(fid, line);

		// get box size;
		for (int i = 0; i != 3; ++i) {
			std::getline(fid, line);
			ss.str(line);
			ss >> Lbox[i][0] >> Lbox[i][1];
			ss.clear();
		}
		std::getline(fid, line);

		// get data;
		for (int i = 0; i != natoms; ++i) {
			std::getline(fid, line);
			ss.str(line);
			for (int j = 0; j != ncolumn; ++j) {
				ss >> data[j][i];
			}
			ss.clear();
		}

#ifdef _WIN32
		printf("\r");
#else
		printf("\r\033[k");
#endif // _WIN32
		fflush(stdout);
		printf("Current Time Step: %d", curTime);

		return true;
	}
	else {
		printf("\n\n### Analysis Done!\n\n"); 
		return false;
	}	
}

double* Trajectory::getData(std::string sign)
{
	auto n = std::find(headline.begin(), headline.end(), sign) - headline.begin();
	return data[n];
}

int Trajectory::getTime()
{
	return curTime;
}

Trajectory::~Trajectory()
{
	fid.close();
	for (int i = 0; i != ncolumn; ++i)
		delete[] data[i];
	delete[] data;
}



