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
		
		/*!
		 * @brief band情報からdosを計算する
		 * @param mc_sim::brillouin_zone バンドの情報
		 * @return std::vector<std::pair<double, double>> dosの疑似curve
		*/
		static std::vector<std::pair<double, double>> dos_tetrahedron(mc_sim::brillouin_zone& bz);
		
		/*!
		 * @brief band情報から(角周波数-domcp)の最大値を計算する
		 * @param std::vector<std::pair<double, double>> dosの疑似curve
		 * @return std::vector<std::pair<double, double>> 温度-domcpmaxの疑似curve
		*/
		static std::vector<std::pair<double, double>> domcpmax_tetrahedron(const std::vector<std::pair<double, double>>&);
		
		/*! シリコンの格子定数 */
		constexpr static inline double si_lattice_constant = 5.431e-10;
	public:
		/*! 比熱などを計算するときに考慮する最大の温度 */
		static inline const int heatcaps_tempmax = 1200;
		
		/*!
		 * @brief 波数空間上の規格座標からLAバンドでの角周波数を<1, 0, 0>モデルで求める関数
		 * @param std::tuple<double, double, double>& 波数空間上の規格座標
		 * @return LAバンド, <1, 0, 0>モデルでの角周波数
		*/
		static double si_angfreq_100_la(const std::tuple<double, double, double>&);
		
		/*!
		 * @brief 波数空間上の規格座標からTAバンドでの角周波数を<1, 0, 0>モデルで求める関数
		 * @param std::tuple<double, double, double>& 波数空間上の規格座標
		 * @return TAバンド, <1, 0, 0>モデルでの角周波数
		*/
		static double si_angfreq_100_ta(const std::tuple<double, double, double>&);
		
		//分散関係テーブルを構築する
		//状態密度テーブル構築他,シミュに前提として要求される(より良い実装求む)
		//static void Si_dispersion_table_construct();
		
		/*!
		 * @brief 分散関係からdosのcurveと, domcpの最大値curveを四面体法で構築する
		 * @param mc_sim::brillouin_zone const & 分散関係
		 * @param std::shared_ptr<mc_sim::logger>& curveに入れるロガー
		 * @return std::pair<curve, curve> dosの関数curveとdomcp最大値の関数curve
		*/
		static std::pair<curve, curve> dos_domcpmax_tetrahedron(mc_sim::brillouin_zone&, std::shared_ptr<mc_sim::logger>&);
		
		/*!
		 * @brief band情報から比熱curveを計算する
		 * @param std::vector<std::shared_ptr<band>> band情報たち
		 * @param std::shared_ptr<mc_sim::logger>& curveに入れるlogger
		 * @return curve 比熱の関数curve
		*/
		static curve heatcap_curve_construct(std::vector<std::shared_ptr<band>>, std::shared_ptr<mc_sim::logger>&);
		
		/*!
		 * @brief band情報から温度とエネルギー体積密度の相互変換を行えるcurveを計算する
		 * @param std::vector<std::shared_ptr<band>> band情報たち
		 * @param std::shared_ptr<mc_sim::logger>& curveに入れるlogger
		 * @return std::pair<curve, curve> 温度からエネルギー体積密度を得るcurveとエネルギー体積密度から温度を得るcurveのpair
		*/
		static std::pair<curve, curve> internal_energy_construct(std::vector<std::shared_ptr<band>>&, std::shared_ptr<mc_sim::logger>&);
		
		/* static double Si_group_velocity(double angular_frequency, int bandnum); */
};
