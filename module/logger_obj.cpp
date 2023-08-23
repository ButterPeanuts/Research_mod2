#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <logger_obj.hpp>

using namespace mc_sim;

//別シンクを使うコンストラクタ
logger_obj::logger_obj(std::string logger_name, std::string file_name){
	auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(file_name);
	stdout_sink->set_level(spdlog::level::info);
	#ifdef MC_SIM_DEBUG
		file_sink->set_level(spdlog::level::debug);
	#else
		file_sink->set_level(spdlog::level::info);
	#endif
	this->logger_sinks = {stdout_sink, file_sink};
	this->set_interface(logger_name);
}

//シンクを流用するコンストラクタ
logger_obj::logger_obj(std::string logger_name, std::vector<spdlog::sink_ptr> sinks){
	this->logger_sinks = sinks;
	this->logger_obj::set_interface(logger_name);
}

//シンクを設定済みの状態からインターフェースを作る
void logger_obj::set_interface(std::string logger_name){
	this->logger_interface = std::make_shared<spdlog::async_logger>(logger_name, this->logger_sinks.begin(), this->logger_sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);
	
	this->debug("This logger is born now!");
}

//thisのシンクを流用して作ったロガーを返す
std::unique_ptr<logger> logger_obj::copy_samesink(std::string logger_name){
	logger_obj samesinklogger(logger_name, this->logger_sinks);
	return std::make_unique<logger_obj>(samesinklogger);
}

void logger_obj::debug(std::string data){
	this->logger_interface->debug(data);
}

void logger_obj::info(std::string data){
	this->logger_interface->info(data);
}

void logger_obj::warn(std::string data){
	this->logger_interface->warn(data);
}

void logger_obj::error(std::string data){
	this->logger_interface->error(data);
}

void logger_obj::critical(std::string data){
	this->logger_interface->critical(data);
}
