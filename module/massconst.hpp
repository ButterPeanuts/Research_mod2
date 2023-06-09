#pragma once
#include <vector>
#include <string>
#include<functional>
#include<algorithm>

#include<curve.hpp>
#include<band.hpp>
#include<logger.hpp>
class massconst {
	private:
		massconst();
		
	public:
		//格子定数 どこに使われてるんだろう
		//static const double Si_lattice_constant;
		//比熱などを計算するときに考慮する最大の温度
		static inline const int heatcaps_tempmax = 600;
		//四面体分割数
		static int Ndiv;
		//整列済みのE_Edgeに対して,k空間体積(鎌倉p17)を計算
		//static double k_volume(std::vector<double> E_Edge, double E);
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
		static curve heatcap_curve_construct(std::vector<band>, mc_sim::logger);
		
		/* static double Si_group_velocity(double angular_frequency, int bandnum); */
};
