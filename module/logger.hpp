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

namespace mc_sim{
	class logger{
		private:
	};
}
