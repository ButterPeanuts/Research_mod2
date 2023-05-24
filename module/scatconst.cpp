#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<scatconst.hpp>
#include<logger.hpp>


scatconst::scatconst(mc_sim::logger& newlogger, std::string filename) : logger(newlogger){
	this->table = std::vector<double>();
	this->logger = newlogger;
	this->file_input(filename);
	this->logger.debug("I'll get " + filename);
}

void scatconst::file_input(std::string filename){
	std::ifstream file(filename);
	if (!file) {
		this->logger.warn(filename + "の読み込みに失敗しました");
		return;
	}
	std::string string_buffer;
	while (std::getline(file, string_buffer)) {
		double part = std::stod(string_buffer);
		
		this->table.push_back(part);
	}
	file.close();
	this->logger.info(filename + "の読み込みに成功しました");
}

double scatconst::a(){return this->table[0];}
double scatconst::b(){return this->table[1];}
double scatconst::chi(){return this->table[2];}
double scatconst::xi(){return this->table[3];}
double scatconst::c(){return this->table[4];}
double scatconst::f(){return this->table[5];}
