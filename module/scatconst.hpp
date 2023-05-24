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
	mc_sim::logger& logger;
	
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
	scatconst(mc_sim::logger& newlogger, std::string filename);
	
	double a();
	double b();
	double chi();
	double xi();
	double c();
	double f();
};
