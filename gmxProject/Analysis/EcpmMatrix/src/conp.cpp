#include "conp.h"
#include <fstream>
#include <chrono>
#include <iomanip>
#include <itp/utility>

Conp::Conp(int gc, char** gv) : argc(gc), argv(gv)
{
}

void Conp::readGro()
{
	// load .gro file;
	std::ifstream file(fnm);
	std::stringstream ss;
	std::string line;
	int allatoms;

	std::getline(file, line);
	std::getline(file, line);

	ss.str(line);
	ss >> allatoms;
	ss.clear();

	x.resize(natoms, 3);
	for (int i = 0; i != natoms; ++i)
	{
		std::getline(file, line);
		ss.str(line.substr(20, 24));
		ss >> x(i, 0) >> x(i, 1) >> x(i, 2);
		ss.clear();
	}

	for (int i = 0; i != allatoms - natoms; ++i)
	{
		std::getline(file, line);
	}

	std::getline(file, line);
	ss.str(line);
	ss >> box[0] >> box[1] >> box[2];
	ss.clear();
	for (auto&& i : box)
	{
		i *= 10;
	}

	x *= 10;
	for (int i = 0; i != natoms; ++i)
	{
		x(i, 2) -= box[2] / 2;
	}
}

void Conp::calc_rReal()
{
	Eigen::ArrayX3d xall(natoms * 9, 3);
	Eigen::ArrayXXd ratio(9, 2);
	ratio << 0, 0, 0, 1, 0, -1, 1, 0, 1, 1, 1, -1, -1, 0, -1, 1, -1, -1;
	for (int i = 0; i != natoms; ++i)
	{
		for (int j = 0; j != ratio.rows(); ++j)
		{
			xall(i + natoms * j, 0) = x(i, 0) + box[0] * ratio(j, 0);
			xall(i + natoms * j, 1) = x(i, 1) + box[1] * ratio(j, 1);
			xall(i + natoms * j, 2) = x(i, 2);
		}
	}

	double r;
	for (int i = 0; i != natoms; ++i)
	{
		for (int j = 0; j != natoms * 9; ++j)
		{
			r = std::pow(x(i, 0) - xall(j, 0), 2) +
				std::pow(x(i, 1) - xall(j, 1), 2) +
				std::pow(x(i, 2) - xall(j, 2), 2);
			if ((r < cutoff * cutoff) && r > 1e-7)
			{
				r = std::sqrt(r);
				rReal(i, j % natoms) = (my_erfc(g_ewald * r) - my_erfc(eta * r / sqrt(2))) / r;
			}
		}
	}
}

void Conp::calc_rSelf()
{
	for (int i = 0; i != natoms; ++i)
	{
		rSelf(i, i) = (std::sqrt(2) * eta - 2 * g_ewald) / std::sqrt(itp::pi);
	}
}

void Conp::calc_rSlab()
{
	for (int i = 0; i != natoms; ++i)
	{
		for (int j = 0; j != natoms; ++j)
		{
			if (b3dc)
			{
				//if system is not slab, no need to add slab correction.
				rSlab(i, j) = 4 * itp::pi * x(i, 2) * x(j, 2) / volume;
			}
		}
	}
}

void Conp::calc_rKspace()
{
	for (int i = 0; i != nthread; ++i)
	{
		thread.emplace_back(&Conp::kspaceThreadFunc, this, i, nthread);
	}
}

void Conp::kspaceThreadFunc(int rank, int numProcess)
{
	double unitk[3];
	unitk[0] = 2.0 * itp::pi / box[0];
	unitk[1] = 2.0 * itp::pi / box[1];
	unitk[2] = 2.0 * itp::pi / box[2];

	int kmax = std::max(kxmax, kymax);
	kmax = std::max(kmax, kzmax);

	double gsqxmx = unitk[0] * unitk[0] * kxmax * kxmax;
	double gsqymx = unitk[1] * unitk[1] * kymax * kymax;
	double gsqzmx = unitk[2] * unitk[2] * kzmax * kzmax;

	double gsqmx = std::max(gsqxmx, gsqymx);
	gsqmx = std::max(gsqmx, gsqzmx);
	gsqmx *= 1.00001;

	int k, l, m;
	double sqk;
	double fc; // Fourier coefficient of the Gaussian function used in the Ewald sum

	double g_ewald_sq_inv = 1.0 / (g_ewald * g_ewald);
	double preu = 8.0 * itp::pi / volume;

	// (k,0,0), (0,l,0), (0,0,m)

	for (m = 1; m <= kmax; m++)
	{
		sqk = (m * unitk[0]) * (m * unitk[0]);
		if (sqk <= gsqmx)
		{
			fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
			for (int ii = rank; ii < natoms; ii += numProcess)
			{
				for (int jj = 0; jj < natoms; ++jj)
				{
					rKspace(ii, jj) += fc * std::cos(m * unitk[0] * (x(ii, 0) - x(jj, 0)));
				}
			}
		}

		sqk = (m * unitk[1]) * (m * unitk[1]);
		if (sqk <= gsqmx)
		{
			fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
			for (int ii = rank; ii < natoms; ii += numProcess)
			{
				for (int jj = 0; jj < natoms; ++jj)
				{
					rKspace(ii, jj) += fc * std::cos(m * unitk[1] * (x(ii, 1) - x(jj, 1)));
				}
			}
		}

		sqk = (m * unitk[2]) * (m * unitk[2]);
		if (sqk <= gsqmx)
		{
			fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
			for (int ii = rank; ii < natoms; ii += numProcess)
			{
				for (int jj = 0; jj < natoms; ++jj)
				{
					rKspace(ii, jj) += fc * std::cos(m * unitk[2] * (x(ii, 2) - x(jj, 2)));
				}
			}
		}
		count += 3;
	}

	// 1 = (k,l,0), 2 = (k,-l,0)

	for (k = 1; k <= kxmax; k++)
	{
		for (l = 1; l <= kymax; l++)
		{
			sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[1] * l) * (unitk[1] * l);
			if (sqk <= gsqmx)
			{
				fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
				for (int ii = rank; ii < natoms; ii += numProcess)
				{
					for (int jj = 0; jj < natoms; ++jj)
					{
						rKspace(ii, jj) += 2 * fc *
							std::cos(k * unitk[0] * (x(ii, 0) - x(jj, 0))) *
							std::cos(l * unitk[1] * (x(ii, 1) - x(jj, 1)));

					}
				}
			}

			count++;
		}
	}

	// 1 = (0,l,m), 2 = (0,l,-m)

	for (l = 1; l <= kymax; l++)
	{
		for (m = 1; m <= kzmax; m++)
		{
			sqk = (unitk[1] * l) * (unitk[1] * l) + (unitk[2] * m) * (unitk[2] * m);
			if (sqk <= gsqmx)
			{
				fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
				for (int ii = rank; ii < natoms; ii += numProcess)
				{
					for (int jj = 0; jj < natoms; ++jj)
					{
						rKspace(ii, jj) += 2 * fc *
							std::cos(l * unitk[1] * (x(ii, 1) - x(jj, 1))) *
							std::cos(m * unitk[2] * (x(ii, 2) - x(jj, 2)));

					}
				}
			}
			count++;
		}
	}

	// 1 = (k,0,m), 2 = (k,0,-m)

	for (k = 1; k <= kxmax; k++)
	{
		for (m = 1; m <= kzmax; m++)
		{
			sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[2] * m) * (unitk[2] * m);
			if (sqk <= gsqmx)
			{
				fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
				for (int ii = rank; ii < natoms; ii += numProcess)
				{
					for (int jj = 0; jj < natoms; ++jj)
					{
						rKspace(ii, jj) += 2 * fc *
							std::cos(k * unitk[0] * (x(ii, 0) - x(jj, 0))) *
							std::cos(m * unitk[2] * (x(ii, 2) - x(jj, 2)));

					}
				}
			}
			count++;
		}
	}

	// 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

	for (k = 1; k <= kxmax; k++)
	{
		for (l = 1; l <= kymax; l++)
		{
			for (m = 1; m <= kzmax; m++)
			{
				sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[1] * l) * (unitk[1] * l) +
					(unitk[2] * m) * (unitk[2] * m);
				if (sqk <= gsqmx)
				{
					fc = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
					for (int ii = rank; ii < natoms; ii += numProcess)
					{
						for (int jj = 0; jj < natoms; ++jj)
						{
							rKspace(ii, jj) += 4 * fc *
								std::cos(k * unitk[0] * (x(ii, 0) - x(jj, 0))) *
								std::cos(l * unitk[1] * (x(ii, 1) - x(jj, 1))) *
								std::cos(m * unitk[2] * (x(ii, 2) - x(jj, 2)));
						}
					}
				}
				count++;
			}
		}
	}
}

void Conp::calc_Matrix()
{
	// get comment line args;
	itp::Getopt getopt(argc, argv, "calculate Matrix for ECPM\n");

	getopt(fnm, "-f", true, "gro file name");
	getopt(natoms, "-n", true, "number of total electrode atoms");
	getopt(kxmax, "-kx", true, "kxmax");
	getopt(kymax, "-ky", true, "kymax");
	getopt(kzmax, "-kz", true, "kzmax");
	getopt(b3dc, "-3dc", true, "3d(0) or 3dc(1)");

	// optional 
	nthread = std::thread::hardware_concurrency();
	getopt(nthread, "-thread", false, "number of thread");
	getopt(cutoff, "-cutoff", false, "cut off(nm)");
	getopt(eta, "-eta", false, "eta");
	getopt(g_ewald, "-ewald", false, "g_ewald");
	getopt.finish();

	// print input args
	fmt::print("\ninput args:\n");
	fmt::print(".gro file name: {}\n", fnm);
	fmt::print("number of total electrode atoms: {}\n", natoms);
	fmt::print("kxmax: {}\n", kxmax);
	fmt::print("kymax: {}\n", kymax);
	fmt::print("kzmax: {}\n", kzmax);
	fmt::print("3d(0) or 3dc(1): {}\n", b3dc);
	fmt::print("number of thread: {}\n", nthread);
	fmt::print("cutoff: {}\n\n", cutoff);

	readGro();

	cutoff *= 10; // transfer nm to A;
	rMatrix.resize(natoms, natoms);
	rReal.resize(natoms, natoms);
	rSelf.resize(natoms, natoms);
	rSlab.resize(natoms, natoms);
	rKspace.resize(natoms, natoms);

	rMatrix.fill(0);
	rReal.fill(0);
	rSelf.fill(0);
	rSlab.fill(0);
	rKspace.fill(0);

	volume = box[0] * box[1] * box[2];

	count = 0;
	calc_rKspace();
	calc_rReal();
	calc_rSelf();
	calc_rSlab();

	// print remained time info:
	std::thread(&Conp::printProgressInfo, this).detach();

	for (auto&& i : thread)
	{
		if (i.joinable())
			i.join();
	}

	rMatrix = rReal + rSlab + rSelf + rKspace;
	rInvMatrix = rMatrix.matrix().inverse().array();

	saveToFile("rMatrix.dat", rMatrix, "{:15.8e} ");
	saveToFile("rInvMatrix.dat", rInvMatrix, "{:15.8e} ");
}

void Conp::get_getPotFile()
{
	// get comment line args;
	itp::Getopt getopt(argc, argv, "get \"getPotFile\" for ECPM\n");

	int allatoms = 0;
	std::vector<double> voltage = { 0, 1, 2 };
	double lowPos = 0.0;
	double upPos = 10.0;

	getopt(fnm, "-f", true, "gro file name");
	getopt(natoms, "-n", true, "number of total electrode atoms");
	getopt(allatoms, "-N", true, "number of total atoms of system");
	getopt.getArray(voltage, "-v", true, "voltage");
	getopt(lowPos, "-low", true, "bulk low position(z)(nm)");
	getopt(upPos, "-up", true, "bulk up position(z)(nm)");

	getopt(cutoff, "-cutoff", false, "cut off(nm)");

	getopt.finish();

	readGro();

	for (auto&& v : voltage)
	{
		std::ofstream file("getPot_parameters_{}V.dat"_format(v));
		fmt::print(file, "NAllatom: {:8d}\n", allatoms);
		fmt::print(file, "rcutoff: {:8.3f}\n", cutoff);
		fmt::print(file, "freqCal: {:8d}\n", 1);
		fmt::print(file, "freqOut: {:8d}\n", 10);
		fmt::print(file, "boundraylow1: {:10.4f}\n", lowPos);
		fmt::print(file, "boundrayup1: {:10.4f}\n", upPos);
		fmt::print(file, "boundraylow2: {:10.4f}\n", lowPos);
		fmt::print(file, "boundrayup2: {:10.4f}\n", upPos);
		fmt::print(file, "NUserSelectGrid: {:8d}\n", natoms);
		for (int i = 0; i != natoms / 2; ++i)
		{
			fmt::print(file, "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:8d}\n", x(i, 0) / 10, x(i, 1) / 10,
				(x(i, 2) + box[2] / 2) / 10, v / 2, i);
		}
		for (int i = natoms / 2; i != natoms; ++i)
		{
			fmt::print(file, "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:8d}\n", x(i, 0) / 10, x(i, 1) / 10,
				(x(i, 2) + box[2] / 2) / 10, -v / 2, i);
		}
		file.close();
	}
}

double Conp::my_erfc(double grij)
{
	constexpr double EWALD_P = 0.3275911;
	constexpr double A1 = 0.254829592;
	constexpr double A2 = -0.284496736;
	constexpr double A3 = 1.421413741;
	constexpr double A4 = -1.453152027;
	constexpr double A5 = 1.061405429;

	double expm2 = std::exp(-grij * grij);
	double t = 1.0 / (1.0 + EWALD_P * grij);
	return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
}

void Conp::printProgressInfo()
{
	int kmax = std::max(kxmax, kymax);
	kmax = std::max(kmax, kzmax);
	int ncount = kmax * 3 + kxmax * kymax + kymax * kzmax + kxmax * kzmax + kxmax * kymax * kzmax;
	ncount *= nthread;

	auto beginTime = std::chrono::system_clock::now();

	while (true)
	{
		auto endTime = beginTime + (std::chrono::system_clock::now() - beginTime) * ncount / int(count);
		auto t_c = std::chrono::system_clock::to_time_t(endTime);

		auto remain = std::chrono::duration_cast<std::chrono::seconds>(endTime - std::chrono::system_clock::now());

		itp::setScrollOutput();
		if (remain.count() <= 240)
		{
			std::cout << "Progress: {:.2f}%, "_format(100.0 * count / ncount)
				<< "remaining wall clock time: " << remain.count() << " s" << std::string(20, ' ');
		}
		else
		{
			std::cout << "Progress: {:.2f}%, "_format(100.0 * count / ncount)
				<< "will finish at " << std::put_time(std::localtime(&t_c), "%F %T");
		}

		if (count >= ncount)
		{
			fmt::print("\nFinished!\n");
			break;
		}

		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
}

void Conp::saveToFile(std::string fnm, const Eigen::ArrayXXd& data, std::string fmt)
{
	std::ofstream ofile(fnm);
	for (int i = 0; i < data.rows(); i++)
	{
		for (int j = 0; j < data.cols(); j++)
		{
			fmt::print(ofile, fmt, data(i, j));
		}
		fmt::print(ofile, "\n");
	}
	ofile.close();
}
