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

double physconst::bedist2(double ang_freq, double temp, double left_const){
	//10^-4から10^4
	double energy_ratio = physconst::dirac / physconst::boltzmann * ang_freq / temp;
	if (energy_ratio > 100){
		while (energy_ratio > 100){
			left_const *= std::exp(-100);
			energy_ratio -= 100;
		}
		return left_const * std::exp(-energy_ratio);
	} else {
		return 1 / (std::exp(energy_ratio) - 1);
	}
}

std::tuple<double, double, double> physconst::indextostd(std::tuple<int, int, int> const & index, int const ndiv){
	return {static_cast<double>(std::get<0>(index)) / static_cast<double>(ndiv), static_cast<double>(std::get<0>(index)) / static_cast<double>(ndiv), static_cast<double>(std::get<0>(index)) / static_cast<double>(ndiv)};
}


double physconst::eukleideia_metrike(std::tuple<double, double, double> const & coor){
	return std::sqrt(std::get<0>(coor) * std::get<0>(coor) + std::get<1>(coor) * std::get<1>(coor) + std::get<2>(coor) * std::get<2>(coor));
}
