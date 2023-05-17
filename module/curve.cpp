#include<vector>
#include<utility>
#include<string>
#include<algorithm>
#include<cmath>
#include<limits>
#include<curve.hpp>
#include<logger.hpp>


curve::curve(mc_sim::logger& newlogger) : logger(newlogger){
	this->table = std::vector<std::pair<double, double>>();
	this->logger = newlogger;
	this->logger.debug("The curve is born now!");
}

curve::curve(mc_sim::logger& newlogger, std::string filename) : curve(newlogger){
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

void curve::append(double x, double y){
	this->table.push_back({x, y});
}
