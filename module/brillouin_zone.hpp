/*!
 * @file brillouin_zone.hpp
 * @brief ブリルアンゾーン上の分散関係を管理する
 * @author ButterPeanuts
 * @data 2023-06-16
*/
#pragma once

#include <modeenum.hpp>
#include <algorithm>
namespace mc_sim{
	/*!
	 * @brief ブリルアンゾーン上の分散関係を管理するクラスのインターフェース
	 * @details ブリルアンゾーン上の(角波数)座標から, 角周波数を計算, 補間して返すなど
	 * \f$\Lambda\f$点が0, \f$X\f$点が\f$N_{div}\f$とするメッシュ座標系を使用する
	 * より正確には, \f$a\f$を格子定数, \f$k_{max} = \frac{2\pi}{aN_{div}}\f$を\f$X\f$点の角波数, \f$(x, y, z)\f$をメッシュ座標とすると, 角波数座標では\f$\frac{k_{max}}{N_{div}}(x, y, z)\f$となる
	*/
	class brillouin_zone{
		public:
		/*!
		 * @brief 角波数空間におけるブリルアンゾーンの分割数\f$N_{div}\f$
		 * @details k点のうち, \f$\Lambda\f$点と\f$X\f$点の分割数
		*/
		virtual int ndiv_getter() = 0;
		
		/*!
		 * @brief メッシュ座標から記録されている角周波数を返す
		 * @param std::tuple<int, int, int> メッシュ座標系でx, y, zの座標を指定するtuple
		 * @return 角周波数
		*/
		virtual double angfreq_index(std::tuple<int, int, int>) = 0;
		
		/*!
		 * @brief 縦波か横波かを取得する関数
		*/
		virtual wave_direction directions_getter() = 0;
		
		/*!
		 * @brief 音響波か光学波か取得する関数
		*/
		virtual wave_mode mode_getter() = 0;
	};
}
