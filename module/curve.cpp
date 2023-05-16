#include<vector>
#include<utility>
#include<string>
#include<curve.hpp>
#include<logger.hpp>


curve::curve(mc_sim::logger& newlogger) : logger(newlogger){
	this->table = std::vector<std::pair<double, double>>();
	this->logger = newlogger;
	this->logger.info("The curve is born now!");
}

curve::curve(mc_sim::logger& newlogger, std::string filename) : logger(newlogger){
	this->table = std::vector<std::pair<double, double>>();
	this->logger = newlogger;
	this->logger.info("The curve is born now!");
	this->file_input(filename);
}
