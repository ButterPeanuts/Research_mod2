/*!
 * @file logger.hpp
 * @brief "ログクラス"のインターフェース
 * @author ButterPeanuts
 * @date 2023-05-17
*/
#pragma once
/**
 * @brief "ログクラス"のインターフェース spdlogへの依存を避けたい
*/

#include <string>
namespace mc_sim{
	template <class self>
	class logger{
		/*!
		 * @brief ロガーのコピー
		 * @details コピーというほどではないが, 同じシンクを持つロガーを作るためのメンバ
		 * @param logger_name ロガーの表示名
		 * @return mc_sim::logger_obj 同じシンクを持つが, thisとは違うロガー
		*/
		virtual self copy_samesink(std::string logger_name) = 0;
		
		virtual void debug(std::string data) = 0;
		virtual void info(std::string data) = 0;
		virtual void warn(std::string data) = 0;
		virtual void error(std::string data) = 0;
		virtual void critical(std::string data) = 0;
	};
}
