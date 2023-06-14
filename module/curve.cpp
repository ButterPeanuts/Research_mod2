#include<vector>
#include<utility>
#include<string>
#include<algorithm>
#include<cmath>
#include<limits>
#include<fstream>
#include<sstream>
#include<curve.hpp>
#include<logger.hpp>


curve::curve(std::shared_ptr<mc_sim::logger>& newlogger) : logger(newlogger){
	this->table = std::vector<std::pair<double, double>>();
	this->logger = newlogger;
	this->logger->debug("The curve is born now!");
}

curve::curve(std::shared_ptr<mc_sim::logger>& newlogger, std::string filename) : curve(newlogger){
	this->file_input(filename);
}

double curve::itpl_getter(double x){
	auto right = std::lower_bound((this->table).begin(), (this->table).end(), std::make_pair(x, 0.0));
	if (right == this->table.begin()){
		//定義域より左側なら
		return right->second;
	}
	if (right == this->table.end()){
		//定義域より右側なら
		return (right - 1)->second;
	}
	//定義域の中なら
	double segmenttime = (x - (right - 1)->first) / (right->first - (right - 1)->first);
	return std::lerp((right - 1)->second, right->second, segmenttime);
}

double curve::max(){
	//負無限大
	auto maxy = this->table.begin();
	for(auto i = this->table.begin(); i < (this->table.end()); i++){
		if (maxy->second < i->second){
			maxy = i;
		}
	}
	return maxy->second;
}

double curve::left_edge(){
	return this->table.begin()->first;
}

double curve::right_edge(){
	return (this->table.end() - 1)->first;
}

void curve::append(double x, double y){
	this->table.push_back({x, y});
	this->logger->debug(std::to_string(x) + ", " + std::to_string(y) + " is added");
}

void curve::file_output(std::string filename){
	std::ofstream file(filename);
	if (!file) {
		this->logger->warn(filename + "の保存に失敗しました");
		return;
	}
	
	for (auto i = this->table.begin(); i < this->table.end(); i++) {
		file << i->first << "," << i->second << std::endl;
	}
	
	file.close();
	this->logger->info(filename + "の保存に成功しました");
}

void curve::file_input(std::string filename){
	std::ifstream file(filename);
	if (!file) {
		this->logger->warn(filename + "の読み込みに失敗しました");
		return;
	}
	std::string string_buffer;
	while (std::getline(file, string_buffer)) {
		std::string separated_buffer;
		std::istringstream separator(string_buffer);
		
		std::getline(separator, separated_buffer, ',');
		double partx = std::stod(separated_buffer);
		std::getline(separator, separated_buffer, ',');
		double party = std::stod(separated_buffer);
		
		this->table.push_back({partx, party});
	}
	file.close();
	this->logger->info(filename + "の読み込みに成功しました");
}
