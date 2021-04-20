#ifndef __CONP_H__
#define __CONP_H__

#include <itp/core>
#include <itp/getopt>

#include <thread>
#include <atomic>

class Conp
{
public:
	Conp() = default;

	Conp(int gc, char** gv);

	void readGro();

	void calc_rReal();

	void calc_rSelf();

	void calc_rSlab();

	void calc_rKspace();

	void calc_Matrix();

	void get_getPotFile();

private:
	void kspaceThreadFunc(int rank, int numProcess);

	double my_erfc(double value);

	void printProgressInfo();

	void saveToFile(std::string fnm, const Eigen::ArrayXXd& data, std::string fmt="{:15.8e} ");

private:
	int argc;
	char** argv;

	std::string fnm = "wall.gro";
	double eta = 1.979;
	double g_ewald = 0.2602844;
	double cutoff = 1.2; // nm;
	int kxmax = 1;
	int kymax = 1;
	int kzmax = 1;
	int b3dc = 1;
	int   natoms = 0;
	int nthread = 16;
	double volume;
	double box[3];

	Eigen::ArrayX3d x;
	Eigen::ArrayXXd rReal, rSelf, rSlab, rKspace, rMatrix, rInvMatrix;

	std::vector<std::thread> thread;
	std::atomic_int64_t count;
};


#endif // !__CONP_H__
