/*!
 * @file logger.hpp
 * @brief ログクラス 外部ライブラリとの緩衝地帯として
 * @author ButterPeanuts
 * @date 2023-05-17
*/
#pragma once
/**
 * @brief ログクラス 外部ライブラリとの緩衝地帯として
 * @details 外部ライブラリ(現在はspdlog)がもし使えなくなったとき
 * 大規模な書き換えが起こるのを防ぐため 最小限の機能を使えるインターフェース?として
 * 依存性の逆転はできない(外部ライブラリをこちらに合わせるわけにはいかないため)が
 * 依存性の吸収ぐらいならできると期待してる
*/

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <string>
namespace mc_sim{
	class logger{
		private:
			/*! spdlogのロガー このクラスが直接操作するところ */
			std::shared_ptr<spdlog::logger> logger_interface;
			
			/*! spdlogのシンク ロガーからロガーを作るときに使う */
			std::shared_ptr<spdlog::sinks::basic_file_sink_mt> logger_sink;
			
			/*!
			 * @brief コンストラクタ
			 * @param logger_name ロガーの表示名
			 * @param sink ロガーのシンク
			*/
			logger(std::string logger_name, std::shared_ptr<spdlog::sinks::basic_file_sink_mt> sink);
			
			/*!
			 * @brief 設定済みのシンクからロガーインターフェースを設定する
			 * @param logger_name ロガーの表示名
			*/
			void set_interface(std::string logger_name);
		public:
			/*!
			 * @brief コンストラクタ
			 * @param logger_name ロガーの表示名
			 * @param file_name ロガーが書き込むファイル
			*/
			logger(std::string logger_name, std::string file_name);
			
			/*!
			 * @brief ロガーのコピー
			 * @details コピーというほどではないが, 同じシンクを持つロガーを作るためのメンバ
			 * @param logger_name ロガーの表示名
			 * @return mc_sim::logger 同じシンクを持つが, thisとは違うロガー
			*/
			mc_sim::logger copy_samesink(std::string logger_name);
			
			void debug(std::string data);
			void info(std::string data);
			void warn(std::string data);
			void error(std::string data);
			void critical(std::string data);
	};
}
