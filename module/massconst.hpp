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
		//バンド番号規則 TA:0,1 LA:2
		static const double Si_scatter_ATA, Si_scatter_ALA, Si_scatter_BTA, Si_scatter_BLA, Si_scatter_chiTA, Si_scatter_chiLA, Si_scatter_xiTA, Si_scatter_xiLA, Si_scatter_C;
		static const double Si_lattice_constant;
		//比熱などを計算するときに考慮する最大の温度
		static inline const int heatcaps_tempmax = 600;
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
		//static void Si_heatcap_table_construct();
		static curve heatcap_curve_construct(std::vector<band>, mc_sim::logger);
		
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
