/*!
 * @file brillouin_zone_funcobj.hpp
 * @brief ブリルアンゾーン上の分散関係を関数から計算する
 * @details brillouin_zoneのインターフェースに従う
 * @author ButterPeanuts
 * @data 2023-06-19
*/
#pragma once

#include <brillouin_zone.hpp>

#include <algorithm>
#include <functional>
namespace mc_sim{
	/*!
	 * @brief ブリルアンゾーン上の分散関係を関数から計算するクラスのインターフェース
	 * @details ブリルアンゾーン上の(角波数)座標から, 角周波数を計算して返すなど
	 * \f$\Lambda\f$点が0, \f$X\f$点が\f$N_{div}\f$とするメッシュ座標系を使用する
	 * より正確には, \f$a\f$を格子定数, \f$k_{max} = \frac{2\pi}{aN_{div}}\f$を\f$X\f$点の角波数, \f$(x, y, z)\f$をメッシュ座標とすると, 角波数座標では\f$\frac{k_{max}}{N_{div}}(x, y, z)\f$となる
	*/
	class brillouin_zone_funcobj : public brillouin_zone{
		private:
		int ndiv;
		wave_direction direction;
		wave_mode mode;
		
		/*!
		 * @brief 波数空間上の規格座標から角周波数を計算する分散関係関数
		 * @details 波数空間上の規格座標とは, \f$\Lambda\f$点を0, (正の)\f$X\f$点を1とする座標
		 * @param std::tuple<double, double, double> 波数空間上の規格座標
		*/
		std::function<double(std::tuple<double, double, double>)>& disparsion_calc;
		
		public:
		/*!
		 * @brief コンストラクタ
		 * @param int 分割数ndiv
		 * @param std::function<double(std::tuple<double, double, double>)> 波数空間上の規格座標から角周波数を計算する分散関係関数disparsion_calc
		 * @param wave_direction 振動方向
		 * @param wave_mode 振動モード
		*/
		brillouin_zone_funcobj(int, std::function<double(std::tuple<double, double, double>)>&, wave_direction, wave_mode);
		
		/*!
		 * @brief 角波数空間におけるブリルアンゾーンの分割数\f$N_{div}\f$
		 * @details k点のうち, \f$\Lambda\f$点と\f$X\f$点の分割数
		*/
		int ndiv_getter();
		
		/*!
		 * @brief メッシュ座標から記録されている角周波数を返す
		 * @param std::tuple<int, int, int> メッシュ座標系でx, y, zの座標を指定するtuple
		 * @return 角周波数
		*/
		double angfreq_index(std::tuple<int, int, int>);
		
		/*!
		 * @brief 縦波か横波かを取得する関数
		*/
		wave_direction directions_getter();
		
		/*!
		 * @brief 音響波か光学波か取得する関数
		*/
		wave_mode mode_getter();
	};
}
