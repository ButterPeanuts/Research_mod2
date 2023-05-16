#include <spdlog/spdlog.h>
#include <logger_obj.hpp>

using namespace mc_sim;

//別シンクを使うコンストラクタ
logger_obj::logger_obj(std::string logger_name, std::string file_name){
	this->logger_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(file_name);
	this->logger_obj::set_interface(logger_name);
}

//シンクを流用するコンストラクタ
logger_obj::logger_obj(std::string logger_name, std::shared_ptr<spdlog::sinks::basic_file_sink_mt> sink){
	this->logger_sink = sink;
	this->logger_obj::set_interface(logger_name);
}

//シンクを設定済みの状態からインターフェースを作る
void logger_obj::set_interface(std::string logger_name){
	this->logger_interface = std::make_shared<spdlog::logger>(logger_name, this->logger_sink);
	#ifdef MC_SIM_DEBUG
		this->logger_interface->set_level(spdlog::level::debug);
	#else
		this->logger_interface->set_level(spdlog::level::info);
	#endif
	this->info("This logger is born now!");
}

//thisのシンクを流用して作ったロガーを返す
logger_obj& logger_obj::copy_samesink(std::string logger_name){
	logger_obj samesinklogger(logger_name, this->logger_sink);
	return samesinklogger;
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
