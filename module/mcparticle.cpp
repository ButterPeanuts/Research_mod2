#include "mcparticle.hpp"

#include "physconst.hpp"
#include<band.hpp>
#include<logger.hpp>

#include<vector>
#include<cmath>
#include<numbers>
#include<iostream>

using namespace mc_sim;

mc_particle::mc_particle(mc_sim::logger& newlogger, double temperature, std::vector<band> bandinj) : logger(newlogger){
	//バンド情報を設定
	this->banddata = bandinj;
	
	//速度方向のベクトルを作る
	this->velocity_pointing = std::vector<double>(mc_particle::dimension, 0);
	
	//角周波数, バンド, 速度方向といった初期状態を決定
	this->inelastic_scattering(temperature);
	
	//変位の初期状態を決定
	this->position = std::vector<double>(mc_particle::dimension, 0);
}

void mc_particle::nextstep(double dt) {
	double velocity = this->band_current->gvelocity_getter(this->angular_frequency);
	
	for (int i = 0; i < mc_particle::dimension; i++) {
		this->position[i] += dt * this->velocity_pointing[i] * velocity;
	}
}

void mc_particle::boundaryscatter_b(double max_x, double max_y, double max_z) {
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

void mc_particle::scatter(double temperature,double dt,double min_structure) {
	//バンド情報
	auto band = this->band_current;
	
	//フォノン相互散乱(ウムクラップ散乱)
	double tui = band->a() * std::pow(this->angular_frequency, band->chi()) * std::pow(temperature, band->xi()) * std::exp(-band->b() / temperature);
	//欠陥散乱
	double tdi = band->c() * pow(this->angular_frequency, 4);
	//境界散乱A
	double tbi = band->gvelocity_getter(this->angular_frequency) * min_structure * band->f();
	
	//非弾性散乱の確率
	double pnes = 1 - std::exp(-dt * tui);
	//弾性散乱の確率
	double pes = 1 - std::exp(-dt * (tdi + tbi));
	if (1 < pnes + pes){
		//確率がおかしいとき(散乱確率が1以上)は警告を発する
		this->logger.warn("Probability of scattering is 1 or more!");
	}
	
	//散乱の決定, 実行
	std::uniform_real_distribution<> randp(0, 1);
	double scattering_factor = randp(physconst::mtrand);
	if (scattering_factor < pnes){
		this->inelastic_scattering(temperature);
	} else if (scattering_factor < (pnes + pes)){
		this->elastic_scattering();
	}
	return;
}

//弾性散乱(速さ変化なし,速度ベクトル方向変化)
void mc_particle::elastic_scattering() {
	//https://qiita.com/aa_debdeb/items/e416ae8a018692fc07eb も参照のこと
	std::uniform_real_distribution<> randcosth(-1, 1);
	double costh = randcosth(physconst::mtrand);
	double phi = randcosth(physconst::mtrand) * std::numbers::pi;
	double sinth = std::sqrt(1 - costh * costh);
	this->velocity_pointing[0] = sinth * std::cos(phi);
	this->velocity_pointing[1] = sinth * std::sin(phi);
	this->velocity_pointing[2] = costh;
	return;
}

void mc_particle::inelastic_scattering(double temperature){
	//https://github.com/ButterPeanuts/Research_mod2/issues/9
	//このissueの通り, DOSテーブル範囲外の積分は無視できるとして組まれている
	
	//バンド種類用distribution
	std::uniform_int_distribution<> randp(0, this->banddata.size() - 1);
	
	//MC粒子分布関数
	std::vector<std::function<double(double)>> mcp_dists;
	std::for_each(this->banddata.begin(), this->banddata.end(), [temperature, &mcp_dists](band& selectedband) -> void{
		std::function<double(double)> distpart = [temperature, &selectedband](double omega) -> double{
			return physconst::bedist2(omega, temperature, selectedband.dos_getter(omega)  / physconst::dirac / omega);
		};
		mcp_dists.push_back(distpart);
	});
	
	while (true) {
		//どのバンド?
		int pr = randp(physconst::mtrand);
		auto selectedband = (this->banddata.begin() + pr);
		auto selecteddist = (mcp_dists.begin() + pr);
		//結果
		auto result = physconst::vonNeumann_rejection(*selecteddist, selectedband->dos_omega_distribution_getter(), selectedband->dos_distribution_getter());
		
		if (result.first) {
			//採用なら...
			this->angular_frequency = result.second;
			this->band_current = selectedband;
			this->logger.debug("My angular frequency is " + std::to_string(this->angular_frequency));
			this->logger.debug("My band is " + std::to_string(this->band_current - this->banddata.begin()));
			break;
		}
	}
	
	//速度方向変更
	this->elastic_scattering();
	
	return;
}
