/*!
 * @file scatconst.hpp
 * @brief 散乱係数を表すデータのクラス
 * @author ButterPeanuts
 * @date 2023-05-10
*/
#pragma once
#include<vector>
#include<string>
#include<logger.hpp>
#include<modeenum.hpp>

/**
 * @brief 散乱係数を表すデータのクラス
*/
class scatconst{
	private:
	/*! ロガー */
	std::shared_ptr<mc_sim::logger> logger;
	
	/*! 係数を表す */
	std::vector<double> table;
	
	/*!
	 * @brief データをファイルから読み込む
	 * @param filename ファイル名を指定する
	*/
	void file_input(std::string filename);
	
	public:
	/*!
	 * @brief コンストラクタ ファイルからのデータを予め読む
	 * @param filename 読み込むファイル名を指定する
	*/
	scatconst(std::shared_ptr<mc_sim::logger>& newlogger, std::string filename);
	
	double a();
	double b();
	double chi();
	double xi();
	double c();
	double f();
	
	/*!
	 * @brief フォノン-フォノン散乱の緩和時間
	 * @param double 角周波数\f$\omega\f$
	 * @param double 温度\f$T\f$
	*/
	double tau_u_inv(double, double);
	/*!
	 * @brief フォノン-欠陥散乱の緩和時間
	 * @param double 角周波数\f$\omega\f$
	*/
	double tau_d_inv(double);
	/*!
	 * @brief フォノン-境界散乱の緩和時間
	 * @param double 群速度\f$\bar(v)\f$
	*/
	double tau_b_inv(double);
};
