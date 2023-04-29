/*!
 * @file physconst.hpp
 * @brief 物理定数などの計算部分
 * @author ButterPeanuts
 * @date 2023-04-26
*/
#pragma once
#include<vector>
#include<random>
#include<cmath>
#include<algorithm>
#include<iostream>
#include<functional>
/*!
 * @brief 物理定数, 他ヘルパー
 * @details 物理定数やヘルパー関数を収めるクラス
*/
class physconst {
	private:
		physconst();
	public:
		/*! ディラック定数 \f$ \hbar = \frac{h}{2\pi} [J\cdot s]\f$*/
		static inline const double dirac = 1.054571817e-34;
		/*! 円周率 */
		static inline const double pi = 3.14159265358979323846;
		/*! ボルツマン定数 \f$ k_B [J/K] \f$*/
		static inline const double boltzmann = 1.380649e-23;
		/*! 一様乱数発生器 */
		static std::mt19937_64 mtrand;
		/*!
		 * @brief ノイマンの棄却法による乱数発生器
		 * @param (*f)(double) 乱数分布, 確率変数を引数とし, 確率密度を返す
		 * @param xi 確率変数の下限
		 * @param xs 確率変数の上限
		 * @param fm 確率密度の上限
		*/
		static double vonNeumann_rejection(double (*f)(double),double xi, double xs, double fm);
};
class massconst {
	private:
		massconst();
	public:
		//バンド番号規則 TA:0,1 LA:2
		static const double Si_scatter_ATA, Si_scatter_ALA, Si_scatter_BTA, Si_scatter_BLA, Si_scatter_chiTA, Si_scatter_chiLA, Si_scatter_xiTA, Si_scatter_xiLA, Si_scatter_C;
		static const double Si_lattice_constant;
		//比熱などを計算するときに考慮する最大の温度
		static const int heatcaps_Tempmax;
		static std::vector<std::vector<std::vector<std::vector<double>>>> Si_dispersion;
		static std::vector<std::vector<double>> Si_DOS_TA;
		static std::vector<std::vector<double>> Si_DOS_LA;
		//比熱
		static std::vector<double> Si_heatcap;
		static int Ndiv;
		//整列済みのE_Edgeに対して,k空間体積(鎌倉p17)を計算
		static double k_volume(std::vector<double> E_Edge, double E);
		static double Si_angular_wavenumber(std::vector<double> Normalized_angular_wavenumber, int bandnum);

		//分散関係テーブルを構築する
		//状態密度テーブル構築他,シミュに前提として要求される(より良い実装求む)
		static void Si_dispersion_table_construct();

		//状態密度テーブルを構築する
		//シミュに前提として要求される(より良い実装求む)
		static void Si_DOS_table_construct();
		static double Si_DOS_table_construct_tetrahedron(double E, int n);

		//比熱テーブルを構築する
		static void Si_heatcap_table_construct();

		static double Si_group_velocity(double angular_frequency, int bandnum);

		static void DOS_table_output(std::string DOS_TA_filename,std::string DOS_LA_filename,std::vector<std::vector<double>> DOS_TA, std::vector<std::vector<double>> DOS_LA);
		static void dispersion_table_output(std::string dispersion_filename, std::vector<std::vector<std::vector<std::vector<double>>>> dispersion);
		static void DOS_table_input(std::string DOS_TA_filename, std::string DOS_LA_filename, std::vector<std::vector<double>>& DOS_TA, std::vector<std::vector<double>>& DOS_LA);
		static void dispersion_table_input(std::string dispersion_filename, std::vector<std::vector<std::vector<std::vector<double>>>> &dispersion);
		static double DOS_interpolation(std::vector<std::vector<double>> DOS, double omega);
		static double Heat_cap_integration(std::function<double(double)> DOS, double a, double b, int Temperature);
		static void heatcap_table_output(std::string heatcap_filename, std::vector<double> heatcap);
		static void heatcap_table_input(std::string heatcap_filename, std::vector<double>& heatcap);
};
