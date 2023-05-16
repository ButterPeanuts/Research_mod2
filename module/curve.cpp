#include<vector>
#include<utility>
#include<string>
#include<curve.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

curve::curve(){
	this->table = std::vector<std::pair<double, double>>();
	auto logger = spdlog::basic_logger_mt("logger", "log/log.log");
	logger -> set_level(spdlog::level::debug);
	std::cout << "Hello, World" << std::endl;
	std::cout << "BUILD_TYPE=" BUILD_TYPE << std::endl;
	spdlog::info("Hello, World");
	logger->info("Hello, World");
	logger->debug("That is japanse true KUSA movie.");
}

curve::curve(std::string filename){
	this->table = std::vector<std::pair<double, double>>();
	this->file_input(filename);
}
