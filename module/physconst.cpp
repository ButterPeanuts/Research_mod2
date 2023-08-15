#include"physconst.hpp"
#include<random>
#include<utility>
#include<functional>
#include<exception>
static std::random_device thrand;

std::mt19937_64 physconst::mtrand(thrand());
std::pair<bool, double> physconst::vonNeumann_rejection(std::function<double(double)>& f, std::uniform_real_distribution<>& xdist, std::uniform_real_distribution<>& fdist, std::mt19937_64& engine){
	double xr = xdist(engine);
	double fr = fdist(engine);
	if (fr <= f(xr))return {true, xr};
	else return {false, 0};
}

double physconst::bedist(double ang_freq, double temp){
	return 1 / (std::exp(physconst::dirac * ang_freq / physconst::boltzmann / temp) - 1);
}

double physconst::bedist2(double ang_freq, double temp, double left_const){
	//10^-4から10^4
	if (temp == 0 || ang_freq == 0)throw std::domain_error("角周波数, または温度が0です");
	double energy_ratio = physconst::dirac / physconst::boltzmann * ang_freq / temp;
	if (energy_ratio > 100){
		while (energy_ratio > 100){
			left_const *= std::exp(-100);
			energy_ratio -= 100;
		}
		return left_const * std::exp(-energy_ratio);
	} else {
		return left_const / (std::exp(energy_ratio) - 1);
	}
}

std::tuple<double, double, double> physconst::indextostd(std::tuple<int, int, int> const & index, int const ndiv){
	return {static_cast<double>(std::get<0>(index)) / static_cast<double>(ndiv), static_cast<double>(std::get<1>(index)) / static_cast<double>(ndiv), static_cast<double>(std::get<2>(index)) / static_cast<double>(ndiv)};
}


double physconst::eukleideia_metrike(std::tuple<double, double, double> const & coor){
	return std::sqrt(std::get<0>(coor) * std::get<0>(coor) + std::get<1>(coor) * std::get<1>(coor) + std::get<2>(coor) * std::get<2>(coor));
}
