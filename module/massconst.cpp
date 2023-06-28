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
#include<sstream>



//整列済みのE_Edgeに対して,k空間体積微分(rath(1973))を計算
double massconst::k_volume(std::array<double, 4> const & omega_edge, double omega){
	if (omega <= omega_edge[0]) {
		return 0;
	}
	if (omega_edge[0] < omega and omega <= omega_edge[1]) {
		double omega_r0 = omega - omega_edge[0];
		double omega_10 = omega_edge[1] - omega_edge[0];
		double omega_20 = omega_edge[2] - omega_edge[0];
		double omega_30 = omega_edge[3] - omega_edge[0];
		return 3 * omega_r0 * omega_r0 / omega_10 / omega_20 / omega_30;
	}
	if (omega_edge[1] < omega and omega <= omega_edge[2]) {
		double omega_1r = omega_edge[1] - omega;
		double omega_10 = omega_edge[1] - omega_edge[0];
		double omega_20 = omega_edge[2] - omega_edge[0];
		double omega_30 = omega_edge[3] - omega_edge[0];
		double omega_21 = omega_edge[2] - omega_edge[1];
		double omega_31 = omega_edge[3] - omega_edge[1];
		double f = omega_1r * omega_1r * (omega_20 + omega_31) / omega_21 / omega_31;
		return (3 * omega_10 - 6 * omega_1r - 3 * f) / omega_20 / omega_30;
	}
	if (omega_edge[2] < omega and omega < omega_edge[3]) {
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

double massconst::si_angfreq_100_la(const std::tuple<double, double, double>& k_std){
	const double norm_k = std::sqrt(std::pow(std::get<0>(k_std), 2) + std::pow(std::get<1>(k_std), 2) + std::pow(std::get<2>(k_std), 2));
	constexpr double k_border_la = 0.524;
	constexpr double tilt1 = 8480.0 * 2.0 * std::numbers::pi / massconst::si_lattice_constant;
	constexpr double tilt2 = 4240.0 * 2.0 * std::numbers::pi / massconst::si_lattice_constant;
	if (norm_k <= k_border_la) {
		return tilt1 * norm_k;
	}
	else if(norm_k <= 1.0) {
		return tilt1 * k_border_la + tilt2 * (norm_k - k_border_la);
	} else {
		constexpr double top = tilt1 * k_border_la + tilt2 * (1.0 - k_border_la);
		return top;
	}
}

double massconst::si_angfreq_100_ta(const std::tuple<double, double, double>& k_std){
	const double norm_k = std::sqrt(std::pow(std::get<0>(k_std), 2) + std::pow(std::get<1>(k_std), 2) + std::pow(std::get<2>(k_std), 2));
	constexpr double k_border_ta = 0.403;
	constexpr double tilt = 5860.0 * 2.0 * std::numbers::pi / massconst::si_lattice_constant;
	if (norm_k <= k_border_ta) {
		return tilt * norm_k;
	}
	else {
		constexpr double top = tilt * k_border_ta;
		return top;
	}
}

//改修完了?
curve massconst::doscurve_tetrahedron(mc_sim::brillouin_zone& bz, std::shared_ptr<mc_sim::logger>& logger){
	//仮のカーブ 外部から書き込み可能にしたいためこの形式
	std::vector<std::pair<double, double>> pscurve;
	//角周波数0を最初に入れておく
	pscurve.emplace_back(0.0, 0.0);
	//あとは初項1, 公比1.001[rad/s]の角周波数をカーブに入れる
	double omega_a = 1.0;
	double omega_r = 1.001;
	double omega_max = std::pow(10.0, 14.0);
	for (double omega = omega_a; omega < omega_max; omega *= omega_r){
		pscurve.emplace_back(omega, 0.0);
	}
	
	//ndivを事前に取得しておく
	int ndiv = bz.ndiv_getter();
	
	//四面体の頂点の各周波数を入れるarray
	std::array<double, 4> omega_edge;
	//並列処理用
	std::vector<future<void>> futures;
	
	//積分処理の呼び出し, 加算を行う
	auto dos_integration = [&pscurve, &omega_edge, &futures]() -> void{
		//edgeをソート
		std::sort(omega_edge.begin(), omega_edge.end());
		//dosの変化が起こる角周波数の領域を探索
		auto start = std::upper_bound(pscurve.begin(), pscurve.end(), std::make_pair(omega_edge[0], std::numeric_limits<double>::infinity()));
		auto stop = std::upper_bound(pscurve.begin(), pscurve.end(), std::make_pair(omega_edge[3], std::numeric_limits<double>::infinity()));
		//領域内の数値に対して, 積分値を求め値を更新するfutureを作る
		std::for_each(start, stop, [&omega_edge, &futures](auto& s){
			futures.push_back(std::async(std::launch::deferred, [&omega_edge, &s](){
				s.second += massconst::k_volume(omega_edge, s.first);
			}));
		});
		//完了を待つ
		for (future<void>& f: futures){
			f.get();
		}
		//futureを消す
		futures.clear();
		return;
	};
	for (int i2 = 0; i2 < ndiv; i2++) {
		for (int j2 = 0; j2 < ndiv; j2++) {
			for (int k2 = 0; k2 < ndiv; k2++) {
				//立方体がブリルアンゾーン外の場合
				if (i2 + j2 + k2 + 3 >= ndiv * 3 / 2 + 3)continue;
				//type1
				{
					omega_edge[0] = bz.angfreq_index({i2, j2, k2});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
					omega_edge[2] = bz.angfreq_index({i2, j2 + 1, k2});
					omega_edge[3] = bz.angfreq_index({i2, j2, k2 + 1});
					dos_integration();
				}
				
				//type1のみブリルアンゾーン内の場合
				if (i2 + j2 + k2 + 3 >= ndiv * 3 / 2 + 2)continue;
				//type2
				{
					omega_edge[0] = bz.angfreq_index({i2, j2, k2 + 1});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
					omega_edge[2] = bz.angfreq_index({i2, j2 + 1, k2 + 1});
					omega_edge[3] = bz.angfreq_index({i2 + 1, j2, k2 + 1});
					dos_integration();
				}
				//type3
				{
					omega_edge[0] = bz.angfreq_index({i2, j2, k2 + 1});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
					omega_edge[2] = bz.angfreq_index({i2, j2 + 1, k2 + 1});
					omega_edge[3] = bz.angfreq_index({i2, j2 + 1, k2});
					dos_integration();
				}
				//type5
				{
					omega_edge[0] = bz.angfreq_index({i2 + 1, j2 + 1, k2});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
					omega_edge[2] = bz.angfreq_index({i2, j2 + 1, k2});
					omega_edge[3] = bz.angfreq_index({i2, j2 + 1, k2 + 1});
					dos_integration();
				}
				//type6
				{
					omega_edge[0] = bz.angfreq_index({i2 + 1, j2 + 1, k2});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2, k2});
					omega_edge[2] = bz.angfreq_index({i2 + 1, j2, k2 + 1});
					omega_edge[3] = bz.angfreq_index({i2, j2 + 1, k2 + 1});
					dos_integration();
				}
				//type4のみブリルアンゾーン外の場合
				if (i2 + j2 + k2 + 3 >= ndiv * 3 / 2 + 1)continue;
				{
					omega_edge[0] = bz.angfreq_index({i2 + 1, j2 + 1, k2});
					omega_edge[1] = bz.angfreq_index({i2 + 1, j2 + 1, k2 + 1});
					omega_edge[2] = bz.angfreq_index({i2 + 1, j2, k2 + 1});
					omega_edge[3] = bz.angfreq_index({i2, j2 + 1, k2 + 1});
					dos_integration();
				}
			}
		}
	}
	
	curve dos(logger);
	const double intconst = 1 / (massconst::si_lattice_constant * massconst::si_lattice_constant * massconst::si_lattice_constant * 6.0 * ndiv * ndiv * ndiv);
	for (auto i: pscurve){
		
		dos.append(i.first, i.second * intconst);
	}
	return dos;
}

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
		if (omega == 0){
			return 0.0;
		};
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

std::pair<curve, curve> massconst::internal_energy_construct(std::vector<std::shared_ptr<band>>& banddata, std::shared_ptr<mc_sim::logger>& newlogger) {
	//内部エネルギーの計算
	curve internal_energy(newlogger);
	curve temperature(newlogger);
	internal_energy.append(0.0, 0.0);
	temperature.append(0.0, 0.0);
	
	//被積分関数
	auto calculator = [](double omega, band& target_band, double t){
		if (t == 0){
			return 0.0;
		}
		if (omega == 0){
			//基本的に0.0
			return target_band.dos_getter(0.0) / physconst::boltzmann / t;
		};
		
		constexpr double diracpark = physconst::dirac / physconst::boltzmann;
		//エネルギー比, これが60を超えたあたりがボルツマン近似域
		//700を超えたあたりがdouble型範囲超過粋
		//0を超えず, 760ぐらいまではある
		double energy_ratio = diracpark * omega / t;
		//被積分関数からボースアインシュタイン統計の部分をぬいたもの
		//大体10 ^ -10ぐらい
		double other = physconst::dirac * omega * target_band.dos_getter(omega);
		if (60.0 < energy_ratio){
			//温度が低いとき, 角周波数が高いとき
			while (60.0 < energy_ratio){
				other *= std::exp(-60.0);
				energy_ratio -= 60.0;
			}
			other *= std::exp(-energy_ratio);
		} else {
			double exp_er = std::exp(energy_ratio);
			other /= (exp_er - 1);
		}
		
		return other;
	};
	std::vector<std::future<std::pair<int, double>>> futures;
	for (int t = 1; t < massconst::heatcaps_tempmax + 1; t++) {
		futures.push_back(std::async(std::launch::async, [t, &banddata, &calculator](){
			//温度tにおける最終的な値を出すlambda
			double energy = 0;
			for (std::shared_ptr<band> i: banddata){
				energy += Romberg(i->dos_leftedge(), i->dos_rightedge(), 10, 10, std::bind(calculator, std::placeholders::_1, std::ref(*i), t));
			}
			return std::make_pair(t, energy);
		}));
	}
	for (auto& i: futures){
		auto res = i.get();
		internal_energy.append(res.first, res.second);
		temperature.append(res.second, res.first);
	}
	return {internal_energy, temperature};
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
