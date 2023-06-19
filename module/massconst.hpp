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
		static double k_volume(std::tuple<double, double, double, double> const &, double omega);
		
	public:
		//格子定数 どこに使われてるんだろう
		//static const double Si_lattice_constant;
		/*! 比熱などを計算するときに考慮する最大の温度 */
		static inline const int heatcaps_tempmax = 1200;
		
		//static double Si_angular_wavenumber(std::vector<double> Normalized_angular_wavenumber, int bandnum);
		
		//分散関係テーブルを構築する
		//状態密度テーブル構築他,シミュに前提として要求される(より良い実装求む)
		//static void Si_dispersion_table_construct();
		
		//状態密度テーブルを構築する
		//シミュに前提として要求される(より良い実装求む)
		//static void Si_DOS_table_construct();
		//static double Si_DOS_table_construct_tetrahedron(double E, int n);
		
		//比熱テーブルを構築する
		//static void Si_heatcap_table_construct();
		static curve heatcap_curve_construct(std::vector<std::shared_ptr<band>>, std::shared_ptr<mc_sim::logger>&);
		
		/* static double Si_group_velocity(double angular_frequency, int bandnum); */
};
