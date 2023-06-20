#include"massconst.hpp"
#include"physconst.hpp"
#include"Integral.h"
#include<curve.hpp>
#include<brillouin_zone.hpp>

#include<vector>
#include<omp.h>
#include<chrono>
#include<functional>
#include<numbers>
#include<future>
#include<memory>
#include<array>
#include<limits>
#include<algorithm>
#include<utility>
#include<future>



//整列済みのE_Edgeに対して,k空間体積微分(rath(1973))を計算
double massconst::k_volume(std::array<double, 4> const & omega_edge, double omega){
	if (omega <= omega_edge[0]) {
		return 0;
	}
	if (omega_edge[0] < omega <= omega_edge[1]) {
		double omega_r0 = omega - omega_edge[0];
		double omega_10 = omega_edge[1] - omega_edge[0];
		double omega_20 = omega_edge[2] - omega_edge[0];
		double omega_30 = omega_edge[3] - omega_edge[0];
		return 3 * omega_r0 * omega_r0 / omega_10 / omega_20 / omega_30;
	}
	if (omega_edge[1] < omega <= omega_edge[2]) {
		double omega_1r = omega_edge[1] - omega;
		double omega_10 = omega_edge[1] - omega_edge[0];
		double omega_20 = omega_edge[2] - omega_edge[0];
		double omega_30 = omega_edge[3] - omega_edge[0];
		double omega_21 = omega_edge[2] - omega_edge[1];
		double omega_31 = omega_edge[3] - omega_edge[1];
		double f = omega_1r * omega_1r * (omega_20 + omega_31) / omega_21 / omega_31;
		return (3 * omega_10 - 6 * omega_1r - 3 * f) / omega_20 / omega_30;
	}
	if (omega_edge[2] < omega < omega_edge[3]) {
		double omega_3r = omega_edge[3] - omega;
		double omega_30 = omega_edge[3] - omega_edge[0];
		double omega_31 = omega_edge[3] - omega_edge[1];
		double omega_32 = omega_edge[3] - omega_edge[2];
		return 3 * omega_3r * omega_3r / omega_30 / omega_31 / omega_32;
	}
	else {
		return 0;
	}
};

/*
double massconst::Si_angular_wavenumber(std::vector<double> Normalized_angular_wavenumber, int bandnum){
	const double Norm_kr = std::sqrt(pow(Normalized_angular_wavenumber[0], 2) + pow(Normalized_angular_wavenumber[1], 2) + pow(Normalized_angular_wavenumber[2], 2));
	const double k_rTA = 0.403;
	const double k_rLA = 0.524;
	if (bandnum == 2) {
		if (Norm_kr < k_rLA) {
			return ((double)8480 * 2 * physconst::pi * Norm_kr / massconst::Si_lattice_constant);
		}
		else {
			return ((double)8480 * 2 * physconst::pi * k_rLA / massconst::Si_lattice_constant) + ((double)4240 * 2 * physconst::pi * (Norm_kr - k_rLA) / massconst::Si_lattice_constant);
		}
	}else {
		if (Norm_kr < k_rTA) {
			return ((double)5860 * 2 * physconst::pi * Norm_kr / massconst::Si_lattice_constant);
		}
		else {
			return ((double)5860 * 2 * physconst::pi * k_rTA / massconst::Si_lattice_constant);
		}
	}
}
*/
/*
void massconst::Si_DOS_table_construct() {
	auto start = std::chrono::system_clock::now();
	for (int n = 0; n <= 2; n += 2){
		for (int i = 0; i <= massconst::Ndiv; i++) {
			for (int j = 0; j <= i; j++) {
				for (int k = 0; k <= j; k++) {
					//ijkで登録する各周波数massconst::Si_dispersion[i][j][k][n]を決定 そこからエネルギーを決定
					if (n == 0) {
						if (std::find_if(massconst::Si_DOS_TA.begin(), massconst::Si_DOS_TA.end(),
						[&](const auto& r) {
							return r[0] == massconst::Si_dispersion[i][j][k][n];
						}) != massconst::Si_DOS_TA.end()) {
							continue;
						}
					}
					else {
						if (std::find_if(massconst::Si_DOS_LA.begin(), massconst::Si_DOS_LA.end(),
						[&](const auto& r) {
							return r[0] == massconst::Si_dispersion[i][j][k][n];
						}) != massconst::Si_DOS_LA.end()) {
							continue;
						}
					}
					double E = massconst::Si_dispersion[i][j][k][n] * physconst::dirac;
					double VS = massconst::Si_DOS_table_construct_tetrahedron(E, n);
					std::vector<double> DOS;
					DOS.push_back(massconst::Si_dispersion[i][j][k][n]);
					DOS.push_back(VS * 16 / (pow(massconst::Ndiv, 3) * pow(massconst::Si_lattice_constant, 3)));
					//std::cout << "(" << n << ")" << DOS[0] << " , " << DOS[1] << std::endl;
					if (n == 0) {
						massconst::Si_DOS_TA.push_back(DOS);
					}
					else {
						massconst::Si_DOS_LA.push_back(DOS);
					}
				}
			}
		}
	}
	std::sort(massconst::Si_DOS_TA.begin(), massconst::Si_DOS_TA.end(), [](const std::vector<double>& alpha, const std::vector<double>& beta) {return alpha[0] < beta[0]; });
	std::sort(massconst::Si_DOS_LA.begin(), massconst::Si_DOS_LA.end(), [](const std::vector<double>& alpha, const std::vector<double>& beta) {return alpha[0] < beta[0]; });
	for (int n = 0; n <= 2; n += 2) {
		double capE = (n == 0 ? 3.0e+13 : 1.0e+14);
		for (double tempE = 1.0e+10; tempE <= capE; (tempE >= 1.0e+13 ? tempE += 1.0e+12 : tempE *= 1.1)) {
			double VS = massconst::Si_DOS_table_construct_tetrahedron(tempE*physconst::dirac, n);
			std::vector<double> DOS;
			DOS.push_back(tempE);
			//鎌倉(2003)の方式 正しくない?
			//DOS.push_back(VS * 16 / (pow(massconst::Ndiv, 3) * pow(massconst::Si_lattice_constant, 3)));
			//rath(1974)の方式 正しいかも
			DOS.push_back(VS / 6 / (pow(massconst::Ndiv, 3) * pow(2 * std::numbers::pi, 3)));
			//std::cout << "(" << n << ")" << DOS[0] << " , " << DOS[1] << std::endl;
			if (n == 0) {
				massconst::Si_DOS_TA.push_back(DOS);
			}
			else {
				massconst::Si_DOS_LA.push_back(DOS);
			}
		}
	}

	std::sort(massconst::Si_DOS_TA.begin(), massconst::Si_DOS_TA.end(), [](const std::vector<double>& alpha, const std::vector<double>& beta) {return alpha[0] < beta[0]; });
	std::sort(massconst::Si_DOS_LA.begin(), massconst::Si_DOS_LA.end(), [](const std::vector<double>& alpha, const std::vector<double>& beta) {return alpha[0] < beta[0]; });

	auto end = std::chrono::system_clock::now();
	auto time = end - start;
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
}
*/

curve massconst::doscurve_tetrahedron(mc_sim::brillouin_zone& bz, std::shared_ptr<mc_sim::logger>& logger){
	std::vector<std::pair<double, double>> pscurve;
	pscurve.emplace_back(0.0, 0.0);
	for (double omega = 1.0; omega < std::pow(10.0, 14.0); omega *= 1.001){
		pscurve.emplace_back(omega, 0.0);
	}
	
	int ndiv = bz.ndiv_getter();

	std::array<double, 4> omega_edge;
	std::vector<future<void>> futures;
	for (int i2 = 0; i2 < ndiv; i2++) {
		for (int j2 = 0; j2 < ndiv; j2++) {
			for (int k2 = 0; k2 < ndiv; k2++) {
				//立方体がブリルアンゾーン外の場合
				if (i2 + j2 + k2 + 3 >= ndiv * 3 / 2 + 3)continue;
				//type1
				omega_edge[0] = bz.angfreq_index({i2, j2, k2});
				omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
				omega_edge[2] = bz.angfreq_index({i2, j2 + 1, k2});
				omega_edge[3] = bz.angfreq_index({i2, j2, k2 + 1});
				//ここから関数化?
				//pscurve, omega_edgeを渡し, futuresをキャプチャーする?
				std::sort(omega_edge.begin(), omega_edge.end());
				auto start = std::upper_bound(pscurve.begin(), pscurve.end(), std::make_pair(omega_edge[0], std::numeric_limits<double>::infinity()));
				auto stop = std::upper_bound(pscurve.begin(), pscurve.end(), std::make_pair(omega_edge[3], std::numeric_limits<double>::infinity()));
				std::for_each(start, stop, [&omega_edge, &futures](auto s){
					futures.push_back(std::async(std::launch::async, [&omega_edge, &s](){
						s.second += massconst::k_volume(omega_edge, s.first);
					}));
				});
				for (future<void>& f: futures){
					f.get();
				}
				futures.clear();
				
				//type1のみブリルアンゾーン内の場合
				if (i2 + j2 + k2 + 3 >= massconst::Ndiv * 3 / 2 + 2)continue;
				//type2
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2 + 1][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
				//type3
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
				//type5
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2 + 1][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2 + 1][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
				//type6
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2 + 1][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2 + 1][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
				//type4のみブリルアンゾーン外の場合
				if (i2 + j2 + k2 + 3 >= massconst::Ndiv * 3 / 2 + 1)continue;
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2 + 1][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2 + 1][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2 + 1][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2 + 1][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
			}
		}
	}
	return curve;
}

/*
void massconst::Si_dispersion_table_construct() {
	massconst::Si_dispersion.assign(massconst::Ndiv + 1, std::vector<std::vector<std::vector<double>>>(massconst::Ndiv + 1, std::vector<std::vector<double>>(massconst::Ndiv + 1, std::vector<double>(3, 0))));
	//int step = 0;
	auto start = std::chrono::system_clock::now();
	for (int i = 0; i <= massconst::Ndiv; i++) {
		for (int j = 0; j <= massconst::Ndiv; j++) {
#pragma omp parallel for
			for (int k = 0; k <= massconst::Ndiv; k++) {
				for (int band = 0; band < 3; band++) {
					std::vector<double> kr = { (double)i / (double)massconst::Ndiv,(double)j / (double)massconst::Ndiv,(double)k / (double)(massconst::Ndiv) };
					massconst::Si_dispersion[i][j][k][band] = Si_angular_wavenumber(kr, band);
					//std::cout << step++ << std::endl;
					//std::cout << i << "," << j << "," << k << "(" << band << ") : " << massconst::Si_dispersion[i][j][k][band] << std::endl;
				}
			}
#pragma omp barrier
		}
	}
	auto end = std::chrono::system_clock::now();
	auto time = end - start;
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << std::endl;
}
*/

//改修完了?
curve massconst::heatcap_curve_construct(std::vector<std::shared_ptr<band>> banddata, std::shared_ptr<mc_sim::logger>& newlogger) {
	//比熱(Heat_cap)の計算
	curve heatcap(newlogger);
	heatcap.append(0.0, 0.0);
	
	//被積分関数
	//hbar / k_b, 7.638 * 10 ^ -12
	double diracpark = physconst::dirac / physconst::boltzmann;
	//hbar ^ 2 / k_b, 8.055 * 10 ^ -46
	double dddpk = diracpark * physconst::dirac;
	auto calculator = [diracpark, dddpk](double omega, band& target_band, double t){
		//エネルギー比, これが60を超えたあたりがボルツマン近似域
		//700を超えたあたりがdouble型範囲超過粋
		//0を超えず, 760ぐらいまではある
		double energy_ratio = diracpark * omega / t;
		//被積分関数からボースアインシュタイン統計の部分をぬいたもの
		//大体10 ^ -10ぐらい
		double other = dddpk * omega * omega * target_band.dos_getter(omega) / t / t;
		if (60.0 < energy_ratio){
			//温度が低いとき, 角周波数が高いとき
			while (60.0 < energy_ratio){
				other *= std::exp(-60.0);
				energy_ratio -= 60.0;
			}
			other *= std::exp(-energy_ratio);
		} else {
			double exp_er = std::exp(energy_ratio);
			other *= exp_er;
			other /= (exp_er - 1);
			other /= (exp_er - 1);
		}
		return other;
	};
	std::vector<std::future<std::pair<int, double>>> futures;
	for (int t = 1; t < massconst::heatcaps_tempmax + 1; t++) {
		futures.push_back(std::async(std::launch::async, [t, &banddata, calculator](){
			//温度tにおける最終的な値を出すlambda
			double cv = 0;
			for (std::shared_ptr<band> i: banddata){
				cv += Romberg(i->dos_leftedge(), i->dos_rightedge(), 7, 7, [calculator, &i, t](double omega){
					return calculator(omega, *i, t);
				});
			}
			return std::make_pair(t, cv);
		}));
	}
	for (auto& i: futures){
		auto res = i.get();
		heatcap.append(res.first, res.second);
	}
	return heatcap;
}

/*
double massconst::Si_group_velocity(double angular_frequency, int bandnum) {
	const double k_rTA = 0.403;
	const double k_rLA = 0.524;
	if (bandnum == 2) {
		if (angular_frequency < ((double)8480 * 2 * physconst::pi * k_rLA / massconst::Si_lattice_constant)) {
			return 8480;
		}
		else {
			return 4240;
		}
	}
	else {
		if (angular_frequency < ((double)5860 * 2 * physconst::pi * k_rTA / massconst::Si_lattice_constant)) {
			return 5860;
		}
		else {
			return 0;
		}
	}
};
*/
