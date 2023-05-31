#include"physconst.hpp"
#include<random>
#include<utility>
#include<functional>
static std::random_device thrand;

std::mt19937_64 physconst::mtrand(thrand());
std::pair<bool, double> physconst::vonNeumann_rejection(std::function<double(double)> f, std::uniform_real_distribution<> xdist, std::uniform_real_distribution<> fdist){
	double xr = xdist(physconst::mtrand);
	double fr = fdist(physconst::mtrand);
	if (fr <= f(xr))return {true, xr};
	else return {false, 0};
}

double physconst::bedist(double ang_freq, double temp){
	return 1 / (std::exp(physconst::dirac * ang_freq / physconst::boltzmann / temp) - 1);
}
