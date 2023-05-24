#include "mcparticles.hpp"

#include "physconst.hpp"
#include<band.hpp>
#include<logger.hpp>
#include<vector>
#include<cmath>
#include<numbers>
#include<iostream>

//多分いらなくなる
#include "massconst.hpp"
using namespace mc_particles;


MCParticles::MCParticles(mc_sim::logger& newlogger, double Energy, double Temperature, std::vector<band> bandinj) : logger(newlogger){
	//バンド情報を設定
	this->banddata = bandinj;
	
	//初期状態u(速度方向)決定
	//https://qiita.com/aa_debdeb/items/e416ae8a018692fc07eb も参照のこと
	this->velocity_pointing = std::vector<double>(MCParticles::dimension, 0);
	//弾性散乱だが, 内部的には速度方向のリセットである
	this->Elastic_scattering();
	
	//初期状態omega,band決定
	//状態密度の情報はMCParticleごとにどの値を保持するか変わりそうなので
	//メンバーを追加しておく必要がある
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	//
	//https://github.com/ButterPeanuts/Research_mod2/issues/9
	//このissueの通り, DOSテーブル範囲外の積分は無視できるとして組まれている
	std::vector<std::uniform_real_distribution<>> randx;
	std::vector<std::uniform_real_distribution<>> randf;
	auto dist_set = [randx, randf](band& x){
		//仮 bandDOSの独立変数
		std::uniform_real_distribution<> temp_randx(0, 0);
		
	};
	for_each(this->banddata.begin(), this->banddata.end(), dist_set);
	//棄却法のために最大値を出しているところ
	//無駄が多い
	//事前にテーブル形式で持っておくべき
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	double maxdis = 0;
	for (auto i = massconst::Si_DOS_LA.begin(); i < massconst::Si_DOS_LA.end(); i++) {
		double P = (*i)[1] * physconst::dirac * (*i)[0] / Energy / (exp(physconst::dirac * (*i)[0] / physconst::boltzmann / Temperature) - 1);
		if (maxdis > P)maxdis = P;
	}
	
	//ただ単純記憶ではダメそう
	//なにせDOS_interpolationがあるので
	//ていうかDOS_interpolationはバンドクラスに移管できそう
	std::uniform_real_distribution<> randf(0, maxdis);
	std::uniform_int_distribution<> randp(0, 2);
	for (;;) {
		double xr = randx(physconst::mtrand);
		double fr = randf(physconst::mtrand);
		int pr = randp(physconst::mtrand);
		auto P = [=](double omega, int p) {
			if (p == 2) {
				return massconst::DOS_interpolation(massconst::Si_DOS_LA, omega) * omega * physconst::dirac / Energy / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
			}
			else {
				return massconst::DOS_interpolation(massconst::Si_DOS_TA, omega) * omega * physconst::dirac / Energy / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
			}
		};
		if (fr <= P(xr, pr)) {
			this->angular_frequency = xr;
			this->bandnum = pr;
			break;
		}
	}
	
	//初期状態r決定
	this->position[0] = 0;
	this->position[1] = 0;
	this->position[2] = 0;
}

void MCParticles::Nextstep(double dt) {
	//これもバンドクラスに統合したほうが良さそう
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	double velocity = massconst::Si_group_velocity(angular_frequency,bandnum);
	
	for (int i = 0; i < MCParticles::dimension; i++) {
		position[i] += dt * velocity_pointing[i] * velocity;
	}
}

void MCParticles::Boundary_Scatter_B(double max_x, double max_y, double max_z) {
	//ここは気になるときに検証and解説付加でいいのか
	if ((this->position)[1] < 0 || max_y < this->position[1]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[1] * this->velocity_pointing[1]);
		velocity_pointing[0] /= sin_oldtheta;
		velocity_pointing[2] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[1] = (this->position[1] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
	}
	if ((this->position)[0] < 0 || max_x < this->position[0]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[0] * this->velocity_pointing[0]);
		velocity_pointing[1] /= sin_oldtheta;
		velocity_pointing[2] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[0] = (this->position[0] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[1] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
	}
	if ((this->position)[2] < 0 || max_z < this->position[2]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[2] * this->velocity_pointing[2]);
		velocity_pointing[0] /= sin_oldtheta;
		velocity_pointing[1] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[2] = (this->position[2] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[1] *= sin_newtheta;
	}
}

void MCParticles::Scatter(double Temperature,double dt,double min_structure) {
	//境界散乱Bが必要
	//フォノンフォノン散乱
	//バンド番号なんとかしないとね
	std::uniform_real_distribution<> randx(0, 1);
	std::uniform_real_distribution<> randcosth(-1, 1);
	
	//フォノン相互
	double Pu;
	if (this->bandnum == 0 || this->bandnum == 1) {
		Pu = massconst::Si_scatter_ATA * pow(this->angular_frequency,massconst::Si_scatter_chiTA) * pow(Temperature,massconst::Si_scatter_xiTA) * std::exp(massconst::Si_scatter_BTA / (-Temperature));
	}
	else {
		Pu = massconst::Si_scatter_ALA * pow(this->angular_frequency,massconst::Si_scatter_chiLA) * pow(Temperature,massconst::Si_scatter_xiLA) * std::exp(massconst::Si_scatter_BLA / (-Temperature));
	}
	if (randx(physconst::mtrand) <= (1 - exp(-dt * Pu))) {
		std::uniform_real_distribution<> randx(std::min((*(massconst::Si_DOS_LA.begin() + 1))[0],(*(massconst::Si_DOS_TA.begin() + 1))[0]), std::max((*(massconst::Si_DOS_LA.end() - 1))[0],(*(massconst::Si_DOS_TA.end() - 1))[0]));
		double maxdis = 0;
		for (auto i = massconst::Si_DOS_LA.begin(); i < massconst::Si_DOS_LA.end(); i++) {
			double P = Pu * (*i)[1] * physconst::dirac * (*i)[0] / (physconst::dirac * this->angular_frequency)/ (exp(physconst::dirac * (*i)[0] / physconst::boltzmann / Temperature) - 1);
			if (maxdis > P)maxdis = P;
		}
		
		std::uniform_real_distribution<> randf(0, maxdis);
		std::uniform_int_distribution<> randp(0, 2);
		for (;;) {
			double xr = randx(physconst::mtrand);
			double fr = randf(physconst::mtrand);
			int pr = randp(physconst::mtrand);
			auto P = [&](double omega, int p) {
				if (p == 2) {
					return Pu * massconst::DOS_interpolation(massconst::Si_DOS_LA, omega) * omega * physconst::dirac / (physconst::dirac * this->angular_frequency) / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
				}
				else {
					return Pu * massconst::DOS_interpolation(massconst::Si_DOS_TA, omega) * omega * physconst::dirac / (physconst::dirac * this->angular_frequency) / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
				}
			};
			if (fr <= P(xr, pr)) {
				this->angular_frequency = xr;
				this->bandnum = pr;
				break;
			}
		}
		return;
	}
	
	//フォノン欠陥
	double Pd = massconst::Si_scatter_C * pow(angular_frequency, 4);
	if (randx(physconst::mtrand) <= (1 - exp(-dt * Pd))) {
		this->Elastic_scattering();
		return;
	}
	
	//フォノン境界A
	//F,Lは後に
	double Pb = massconst::Si_group_velocity(angular_frequency, bandnum) * min_structure * 0.55;
	if (randx(physconst::mtrand) <= (1 - exp(-dt * Pb))) {
		this->Elastic_scattering();
		return;
	}
}

//弾性散乱(速さ変化なし,速度ベクトル方向変化)
void MCParticles::Elastic_scattering() {
	std::uniform_real_distribution<> randcosth(-1, 1);
	double costh = randcosth(physconst::mtrand);
	double phi = randcosth(physconst::mtrand) * std::numbers::pi;
	double sinth = std::sqrt(1 - costh * costh);
	this->velocity_pointing[0] = sinth * std::cos(phi);
	this->velocity_pointing[1] = sinth * std::sin(phi);
	this->velocity_pointing[2] = costh;
	return;
}
