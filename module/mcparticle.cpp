#include "mcparticle.hpp"

#include "physconst.hpp"
#include<band.hpp>
#include<logger.hpp>

#include<vector>
#include<cmath>
#include<numbers>
#include<iostream>

using namespace mc_sim;

mc_particle::mc_particle(const std::shared_ptr<mc_sim::logger>& newlogger, double temperature, const std::vector<std::shared_ptr<band>>& bandinj, uint_fast64_t seed) : banddata(bandinj), logger(newlogger){
	//乱数器を生成
	this->mtrand = std::mt19937_64(seed);
	
	//速度方向のベクトルを作る
	this->velocity_pointing = std::vector<double>(mc_particle::dimension, 0);
	
	//角周波数, バンド, 速度方向といった初期状態を決定
	this->angfreq_replace(temperature, false);
	this->elastic_scattering();
	
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
		double cos_newtheta = std::sqrt(randR(this->mtrand));
		velocity_pointing[1] = (this->position[1] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
		this->logger->info("Reflec");
	}
	if ((this->position)[0] < 0 || max_x < this->position[0]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[0] * this->velocity_pointing[0]);
		velocity_pointing[1] /= sin_oldtheta;
		velocity_pointing[2] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(this->mtrand));
		velocity_pointing[0] = (this->position[0] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[1] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
		this->logger->info("Reflec");
	}
	if ((this->position)[2] < 0 || max_z < this->position[2]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[2] * this->velocity_pointing[2]);
		velocity_pointing[0] /= sin_oldtheta;
		velocity_pointing[1] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(this->mtrand));
		velocity_pointing[2] = (this->position[2] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[1] *= sin_newtheta;
		this->logger->info("Reflec");
	}
}

void mc_particle::scatter(double temperature,double dt,double min_structure) {
	//バンド情報
	auto band = this->band_current;
	
	/* //フォノン相互散乱(ウムクラップ散乱) */
	/* double tui = band->a() * std::pow(this->angular_frequency, band->chi()) * std::pow(temperature, band->xi()) * std::exp(-band->b() / temperature); */
	/* //欠陥散乱 */
	/* double tdi = band->c() * pow(this->angular_frequency, 4); */
	/* //境界散乱A */
	/* double tbi = band->gvelocity_getter(this->angular_frequency) / min_structure / band->f(); */
	
	//非弾性散乱の確率
	/* double pnes = 1 - std::exp(-dt * tui); */
	double pnes = band->scatprob_u(this->angular_frequency, temperature, dt);
	//弾性散乱の確率
	/* double pes = (1 - std::exp(-dt * tdi)); */
	double pes = band->scatprob_d(this->angular_frequency, dt) + band->scatprob_b(band->gvelocity_getter(this->angular_frequency), min_structure, dt);
	if (1 < pnes + pes){
		//確率がおかしいとき(散乱確率が1以上)は警告を発する
		this->logger->warn("Probability of scattering is 1 or more!");
	}
	
	//散乱の決定, 実行
	std::uniform_real_distribution<> randp(0, 1);
	double scattering_factor = randp(this->mtrand);
	if (scattering_factor < pnes){
		this->inelastic_scattering(temperature, dt);
	} else if (scattering_factor < (pnes + pes)){
		this->elastic_scattering();
	}
	return;
}

//弾性散乱(速さ変化なし,速度ベクトル方向変化)
void mc_particle::elastic_scattering() {
	//https://qiita.com/aa_debdeb/items/e416ae8a018692fc07eb も参照のこと
	std::uniform_real_distribution<> randcosth(-1, 1);
	double costh = randcosth(this->mtrand);
	double phi = randcosth(this->mtrand) * std::numbers::pi;
	double sinth = std::sqrt(1 - costh * costh);
	this->velocity_pointing[0] = sinth * std::cos(phi);
	this->velocity_pointing[1] = sinth * std::sin(phi);
	this->velocity_pointing[2] = costh;
	return;
}

void mc_particle::inelastic_scattering(double temperature, double dt){
	this->angfreq_replace(temperature, true, dt);
	this->elastic_scattering();
}

void mc_particle::angfreq_replace(double temperature, bool kirchhoff, double dt){
	//https://github.com/ButterPeanuts/Research_mod2/issues/9
	//このissueの通り, DOSテーブル範囲外の積分は無視できるとして組まれている
	
	//バンド種類用distribution
	std::uniform_int_distribution<> randp(0, this->banddata.size() - 1);
	
	//MC粒子分布関数
	std::vector<std::function<double(double)>> mcp_dists;
	for (std::shared_ptr<band> selectedband: this->banddata){
		if (kirchhoff){
			mcp_dists.push_back([temperature, selectedband, dt](double omega) -> double{
				return selectedband->scatprob_u(omega, temperature, dt) * physconst::bedist2(omega, temperature, selectedband->dos_getter(omega) * physconst::dirac * omega);
			});
		}else{
			mcp_dists.push_back([temperature, selectedband](double omega) -> double{
				return physconst::bedist2(omega, temperature, selectedband->dos_getter(omega) * physconst::dirac * omega);
			});
		}
	}
	
	auto domcp_distribution = this->max_dd(temperature);
	while (true) {
		//どのバンド?
		int pr = randp(this->mtrand);
		auto selectedband = *(this->banddata.begin() + pr);
		auto selecteddist = (mcp_dists.begin() + pr);
		//結果
		auto result = physconst::vonNeumann_rejection(*selecteddist, selectedband->dos_omega_distribution_getter(), domcp_distribution, this->mtrand);
		
		if (result.first) {
			//採用なら...
			this->angular_frequency = result.second;
			this->band_current = selectedband;
			this->logger->debug("My angular frequency is " + std::to_string(this->angular_frequency));
			this->logger->debug("My band is " + std::to_string(pr));
			break;
		}
	}
	
	return;
}

std::uniform_real_distribution<double> mc_particle::max_dd(double t){
	return (*std::max_element(banddata.begin(), banddata.end(), [t](const std::shared_ptr<band>& a, const std::shared_ptr<band>& b){return (a->domcp_distribution_getter(t).max() <= b->domcp_distribution_getter(t).max());}))->domcp_distribution_getter(t);
}
