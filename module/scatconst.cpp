#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<scatconst.hpp>
#include<logger.hpp>


scatconst::scatconst(std::shared_ptr<mc_sim::logger>& newlogger, std::string filename) : logger(newlogger){
	this->table = std::vector<double>();
	this->logger = newlogger;
	this->file_input(filename);
	this->logger->debug("I'll get " + filename);
}

void scatconst::file_input(std::string filename){
	std::ifstream file(filename);
	if (!file) {
		this->logger->warn(filename + "の読み込みに失敗しました");
		return;
	}
	std::string string_buffer;
	while (std::getline(file, string_buffer)) {
		double part = std::stod(string_buffer);
		
		this->table.push_back(part);
	}
	file.close();
	this->logger->info(filename + "の読み込みに成功しました");
}

double scatconst::a(){return this->table[0];}
double scatconst::b(){return this->table[1];}
double scatconst::chi(){return this->table[2];}
double scatconst::xi(){return this->table[3];}
double scatconst::c(){return this->table[4];}
double scatconst::f(){return this->table[5];}

double scatconst::tau_u_inv(double omega, double t){
	return this->a() * std::pow(omega, this->chi()) * std::pow(t, this->xi()) * std::exp(-this->b() / t);
}

double scatconst::tau_d_inv(double omega){
	return this->c() * std::pow(omega, 4);
}

double scatconst::tau_b_inv(double gvelocity, double l){
	return gvelocity / this->f() / l;
}
