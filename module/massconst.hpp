/*!
 * @file massconst.hpp
 * @brief 物理定数, 物理関数計算に関するヘルパー
 * @author ButterPeanuts
 * @data 2023-06-19
*/
#pragma once
#include <vector>
#include <string>
#include<functional>
#include<algorithm>
#include<memory>

#include<curve.hpp>
#include<band.hpp>
#include<logger.hpp>
#include<brillouin_zone.hpp>

/*!
 * @brief 物理定数, 物理関数計算に関するヘルパー
*/
class massconst {
	private:
		massconst();
		
		/*!
		 * @brief 整列済みの四面体頂点における角周波数に対して, k空間体積微分(rath(1973))を計算
		 * @param std::tuple<double, double, double, dobule> 四面体頂点における角周波数
		 * @param double 求めたい面の角周波数
		 * @return 四面体内規格化体積
		*/
		static double k_volume(std::array<double, 4>const &, double omega);
		
	public:
		/*! 比熱などを計算するときに考慮する最大の温度 */
		static inline const int heatcaps_tempmax = 1200;
		
		//static double Si_angular_wavenumber(std::vector<double> Normalized_angular_wavenumber, int bandnum);
		
		//分散関係テーブルを構築する
		//状態密度テーブル構築他,シミュに前提として要求される(より良い実装求む)
		//static void Si_dispersion_table_construct();
		
		/*!
		 * @brief 分散関係からdosのcurveを四面体法で構築する
		 * @param mc_sim::brillouin_zone const & 分散関係
		 * @param std::shared_ptr<mc_sim::logger>& curveに入れるロガー
		 * @return curve dosの関数curve
		*/
		static curve doscurve_tetrahedron(mc_sim::brillouin_zone&, std::shared_ptr<mc_sim::logger>&);
		
		/*!
		 * @brief band情報から比熱curveを計算する
		 * @param std::vector<std::shared_ptr<band>> band情報たち
		 * @param std::shared_ptr<mc_sim::logger>& curveに入れるlogger
		 * @return curve 比熱の関数curve
		*/
		static curve heatcap_curve_construct(std::vector<std::shared_ptr<band>>, std::shared_ptr<mc_sim::logger>&);
		
		/* static double Si_group_velocity(double angular_frequency, int bandnum); */
};
