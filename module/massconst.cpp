#include"massconst.hpp"
#include"physconst.hpp"
#include"Integral.h"
#include<curve.hpp>

#include<vector>
#include<omp.h>
#include<chrono>
#include<fstream>
#include<sstream>
#include<string>
#include<functional>
#include<numbers>
#include<future>

const double massconst::Si_scatter_ATA = 1.05e-12;
const double massconst::Si_scatter_ALA = 4.5e-21;
const double massconst::Si_scatter_BTA = 0;
const double massconst::Si_scatter_BLA = 0;
const double massconst::Si_scatter_chiTA = 1;
const double massconst::Si_scatter_chiLA = 2;
const double massconst::Si_scatter_xiTA = 3.86;
const double massconst::Si_scatter_xiLA = 1.36;
const double massconst::Si_scatter_C = 4.95e-45;

const double massconst::Si_lattice_constant = 5.4301e-10;


std::vector<std::vector<std::vector<std::vector<double>>>> massconst::Si_dispersion;
std::vector<std::vector<double>> massconst::Si_DOS_TA;
std::vector<std::vector<double>> massconst::Si_DOS_LA;
std::vector<double> massconst::Si_heatcap;
int massconst::Ndiv;

//整列済みのE_Edgeに対して,k空間体積微分(rath(1973))を計算
double massconst::k_volume(std::vector<double> E_Edge, double E) {
	double E10 = E_Edge[1] - E_Edge[0];
	double E20 = E_Edge[2] - E_Edge[0];
	double E30 = E_Edge[3] - E_Edge[0];
	double E21 = E_Edge[2] - E_Edge[1];
	double E31 = E_Edge[3] - E_Edge[1];
	double E32 = E_Edge[3] - E_Edge[2];
	double f0 = (E - E_Edge[0]) * (E - E_Edge[0]) / (E10 * E20 * E30) / 2;
	double f1 = (E - E_Edge[1]) * (E - E_Edge[1]) / (E10 * E21 * E31) / 2;
	double f3 = (E - E_Edge[3]) * (E - E_Edge[3]) / (E30 * E31 * E32) / 2;
	if (E <= E_Edge[0]) {
		return 0;
	}
	if (E < E_Edge[1]) {
		//return ((E - E_Edge[0]) * (E - E_Edge[0]) / (2 * (E_Edge[1] - E_Edge[0]) * (E_Edge[2] - E_Edge[0]) * (E_Edge[3] - E_Edge[0])));
		return f0;
	}
	if (E <= E_Edge[2]) {
		/*
		double E0 = (E_Edge[2] + E_Edge[3]) * (E_Edge[1] * E_Edge[0] - E * E);
		double E1 = (E_Edge[0] + E_Edge[1]) * (E * E - E_Edge[2] * E_Edge[3]);
		double E2 = 2 * (E_Edge[2] * E_Edge[3] - E_Edge[0] * E_Edge[1]) * E;
		*/
		//return (E0 + E1 + E2) / (2 * E20 * E30 * E21 * E31);
		if (E_Edge[1] == E_Edge[2]) {
			return 1 / (E_Edge[3] - E_Edge[0]) / 6;
		}
		if (E_Edge[0] == E_Edge[1]) {
			double f0u = 2 * (E - E_Edge[1]) * (E_Edge[3] - E) - (E - E_Edge[1]) * (E - E_Edge[1]);
			double f0d = E21 * E31 * E30;
			double f1u = 2 * (E_Edge[2] - E) * (E - E_Edge[0]);
			double f1d = E20 * E21 * E30;
			return (f0u / f0d + f1u / f1d) / 6;
		}
		return f0 - f1;
	}
	if (E < E_Edge[3]) {
		//return (E_Edge[3] - E) * (E_Edge[3] - E) / (2 * E30 * E31 * E32);
		return f3;
	}
	else {
		return 0;
	}
};

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

void massconst::Si_DOS_table_construct() {
	auto start = std::chrono::system_clock::now();
	/*
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
	*/
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

double massconst::Si_DOS_table_construct_tetrahedron(double E, int n) {
	double VS = 0;
#pragma omp parallel for
	for (int i2 = 0; i2 < massconst::Ndiv; i2++) {
		for (int j2 = 0; j2 < massconst::Ndiv; j2++) {
			for (int k2 = 0; k2 < massconst::Ndiv; k2++) {
				std::vector<double> E_Edge;
				//立方体がブリルアンゾーン外の場合
				if (i2 + j2 + k2 + 3 >= massconst::Ndiv * 3 / 2 + 3)continue;
				//type1
				E_Edge.clear();
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2 + 1][j2][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2 + 1][k2][n]);
				E_Edge.push_back(physconst::dirac * massconst::Si_dispersion[i2][j2][k2 + 1][n]);
				std::sort(E_Edge.begin(), E_Edge.end());
				VS += massconst::k_volume(E_Edge, E);
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
#pragma omp barrier
	return VS;
}

void massconst::Si_dispersion_table_construct() {
	massconst::Si_dispersion.assign(massconst::Ndiv + 1, std::vector<std::vector<std::vector<double>>>(massconst::Ndiv + 1, std::vector<std::vector<double>>(massconst::Ndiv + 1, std::vector<double>(3, 0))));
	/* int step = 0; */
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

//改修完了?
curve massconst::heatcap_curve_construct(std::vector<band> banddata, mc_sim::logger newlogger) {
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
	std::vector<std::future<double>> futures;
	for (int t = 1; t < massconst::heatcaps_tempmax + 1; t++) {
		futures.push_back(std::async(std::launch::async, [t, &banddata, calculator](){
			//温度tにおける最終的な値を出すlambda
			double cv = 0;
			for (band& i: banddata){
				cv += Romberg(i.dos_leftedge(), i.dos_rightedge(), 7, 7, [calculator, &i, t](double omega){
					return calculator(omega, i, t);
				});
			}
			return cv;
		}));
	}
	for (auto& i: futures){
		i.get();
	}
	return heatcap;
}

double massconst::Heat_cap_integration(std::function<double(double)> DOS, double a, double b,int Temperature) {
	return Romberg(a, b, 10, 10, [DOS,Temperature](double omega) {
		//当初は久木田(2014)の式をそのまま運用していた
		//Energy_ratioがdouble型の指数上限に引っかかることが判明
		//Energy_ratio^-1 = inv_Energy_ratioを新たに使用する
		/*
		double Energy_ratio = exp(physconst::dirac * omega / physconst::boltzmann / (double)Temperature);
		std::cout << Energy_ratio << std::endl;
		*/
		double inv_Energy_ratio = exp(-1 * physconst::dirac * omega / physconst::boltzmann / (double)Temperature);
		double res = DOS(omega) * pow(physconst::dirac * omega / (double)Temperature, 2) / physconst::boltzmann * inv_Energy_ratio / pow(1 - inv_Energy_ratio, 2);
		//std::cout << res << std::endl;
		return res;
	});
}

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

void massconst::DOS_table_output(std::string DOS_TA_filename, std::string DOS_LA_filename, std::vector<std::vector<double>> DOS_TA, std::vector<std::vector<double>> DOS_LA) {
	std::ofstream DOS_TA_file(DOS_TA_filename + ".txt");
	std::ofstream DOS_LA_file(DOS_LA_filename + ".txt");
	if (!DOS_TA_file || !DOS_LA_file) {
		std::cout << "保存に失敗しました" << std::endl;
		return;
	}
	for (auto i = DOS_TA.begin(); i < DOS_TA.end(); i++) {
		DOS_TA_file << (*i)[0] << "," << (*i)[1] << std::endl;
	}
	for (auto i = DOS_LA.begin(); i < DOS_LA.end(); i++) {
		DOS_LA_file << (*i)[0] << "," << (*i)[1] << std::endl;
	}
	DOS_TA_file.close();
	DOS_LA_file.close();
}

void massconst::dispersion_table_output(std::string dispersion_filename, std::vector<std::vector<std::vector<std::vector<double>>>> dispersion) {

	for (int band = 0; band < 3; band++) {
		std::ofstream dispersion_file(dispersion_filename + std::to_string(band) + ".txt");
		if (!dispersion_file) {
			std::cout << "保存に失敗しました" << std::endl;
			return;
		}
		for (auto i = dispersion.begin(); i < dispersion.end(); i++) {
			dispersion_file << "x[" << i - dispersion.begin() << "]" << std::endl;
			for (auto j = (*i).begin(); j < (*i).end(); j++) {
				for (auto k = (*j).begin(); k < (*j).end(); k++) {
					dispersion_file << (*k)[band];
					if (k < ((*j).end() - 1)) {
						dispersion_file << ",";
					}
					else {
						dispersion_file << std::endl;
					}
				}
			}
		}
	}
}

void massconst::DOS_table_input(std::string DOS_TA_filename, std::string DOS_LA_filename, std::vector<std::vector<double>>& DOS_TA, std::vector<std::vector<double>>& DOS_LA) {
	std::ifstream DOS_TA_file(DOS_TA_filename + ".txt", std::ios::in);
	std::ifstream DOS_LA_file(DOS_LA_filename + ".txt", std::ios::in);
	if (!DOS_TA_file || !DOS_LA_file) {
		std::cout << "保存に失敗しました" << std::endl;
		return;
	}
	std::string string_buffer;
	while (std::getline(DOS_TA_file, string_buffer)) {
		std::string separated_buffer;
		std::vector<double> DOS_buffer;

		std::istringstream separating(string_buffer);

		std::getline(separating, separated_buffer, ',');
		DOS_buffer.push_back(std::stod(separated_buffer));

		std::getline(separating, separated_buffer, ',');
		DOS_buffer.push_back(std::stod(separated_buffer));

		DOS_TA.push_back(DOS_buffer);
	}
	while (std::getline(DOS_LA_file, string_buffer)) {
		std::string separated_buffer;
		std::vector<double> DOS_buffer;

		std::istringstream separating(string_buffer);

		std::getline(separating, separated_buffer, ',');
		DOS_buffer.push_back(std::stod(separated_buffer));

		std::getline(separating, separated_buffer, ',');
		DOS_buffer.push_back(std::stod(separated_buffer));

		DOS_LA.push_back(DOS_buffer);
	}
	DOS_TA_file.close();
	DOS_LA_file.close();
}

void massconst::dispersion_table_input(std::string dispersion_filename, std::vector<std::vector<std::vector<std::vector<double>>>> &dispersion) {
	std::vector<std::ifstream> dispersion_file;
	const int bandnum = 3;
	std::vector<std::string> line_buffer = std::vector<std::string>(bandnum);
	for (int band = 0; band < bandnum; band++) {
		dispersion_file.push_back(std::ifstream(dispersion_filename + std::to_string(band) + ".txt", std::ios::in));
		if (!dispersion_file[band]) {
			std::cout << "保存に失敗しました" << std::endl;
			return;
		}
	}
	std::vector<std::vector<std::vector<double>>> xbuf;
	std::vector<std::vector<double>> ybuf;
	std::vector<double> zbuf;
	bool startskip = true;
	while (1) {
		bool eof = false;
		bool xc = false;
		for (int band = 0; band < bandnum; band++) {
			std::getline(dispersion_file[band], line_buffer[band]);
			if (!dispersion_file[band] || (line_buffer[band][0] == 'x')) {
				xc = true;
			}
			if (!dispersion_file[band])eof = true;
		}
		if (startskip) {
			startskip = false;
			continue;
		}
		if (xc) {
			dispersion.push_back(xbuf);
			xbuf.clear();
		}
		if (eof)break;
		if (xc)continue;
		std::vector<std::istringstream> separating;
		for (int band = 0; band < bandnum; band++)separating.push_back(std::istringstream(line_buffer[band]));
		while (1) {
			bool eof = false;
			for (int band = 0; band < bandnum; band++) {
				std::string separated;
				std::getline(separating[band], separated,',');
				if (!separating[band]) {
					eof = true;
					break;
				}
				zbuf.push_back(stod(separated));
			}
			if (eof)break;
			ybuf.push_back(zbuf);
			zbuf.clear();
		}
		xbuf.push_back(ybuf);
		ybuf.clear();
	}
	for (auto i = dispersion_file.begin(); i < dispersion_file.end(); i++) {
		(*i).close();
	}
}

double massconst::DOS_interpolation(std::vector<std::vector<double>> DOS, double omega) {
	auto match = std::find_if(DOS.begin(), DOS.end(),
	[&](const auto& r) {
		return (r[0] > omega);
	});
	if (match == DOS.begin()) {
		return (*(DOS.begin()))[1];
	}
	if (match == DOS.end()) {
		return (*(DOS.end() - 1))[1];
	}
	double standard = ((*match)[0] - omega) / ((*match)[0] - (*(match - 1))[0]);
	return standard * (*(match - 1))[1] + (1.0 - standard) * (*match)[1];
}

//比熱を保存
void massconst::heatcap_table_output(std::string heatcap_filename ,std::vector<double> heatcap) {
	std::ofstream heatcap_file(heatcap_filename + ".txt");
	if (!heatcap_file) {
		std::cout << "保存に失敗しました" << std::endl;
		return;
	}
	for (auto i = heatcap.begin(); i < heatcap.end(); i++) {
		heatcap_file << (*i) <<std::endl;
	}
	heatcap_file.close();
}

void massconst::heatcap_table_input(std::string heatcap_filename, std::vector<double>& heatcap) {
	std::ifstream heatcap_file(heatcap_filename + ".txt", std::ios::in);
	if (!heatcap_file) {
		std::cout << "読み込みに失敗しました" << std::endl;
		return;
	}
	std::string string_buffer;
	while (std::getline(heatcap_file, string_buffer)) {
		heatcap.push_back(std::stod(string_buffer));
	}
	heatcap_file.close();
}
