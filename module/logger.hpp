/*!
 * @file logger.hpp
 * @brief "ログクラス"のインターフェース
 * @author ButterPeanuts
 * @date 2023-05-17
*/
#pragma once

#include <string>
#include <memory>
namespace mc_sim{
	/**
	 * @brief "ログクラス"のインターフェース spdlogへの依存を避けたい
	*/
	class logger{
		public:
		/*!
		 * @brief ロガーのコピー
		 * @details コピーというほどではないが, 同じシンクを持つロガーを作るためのメンバ
		 * @param logger_name ロガーの表示名
		 * @return mc_sim::logger_obj 同じシンクを持つが, thisとは違うロガー
		*/
		virtual std::unique_ptr<mc_sim::logger> copy_samesink(std::string logger_name) = 0;
		
		virtual void debug(std::string data) = 0;
		virtual void info(std::string data) = 0;
		virtual void warn(std::string data) = 0;
		virtual void error(std::string data) = 0;
		virtual void critical(std::string data) = 0;
	};
}
