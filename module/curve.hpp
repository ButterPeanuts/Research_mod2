/*!
 * @file curve.hpp
 * @brief 二次元グラフ上にプロットできるデータを示すオブジェクトのクラス
 * @author ButterPeanuts
 * @date 2023-05-10
*/
#pragma once
#include<vector>
#include<utility>
#include<string>
#include<logger.hpp>

/**
 * @brief 二次元グラフ上にプロットできるデータを示すオブジェクトのクラス
*/
class curve{
	private:
	/*! プロットデータを表す */
	std::vector<std::pair<double, double>> table;
	
	mc_sim::logger& logger;
	
	public:
	/*!
	 * @brief コンストラクタ
	*/
	curve(mc_sim::logger& newlogger);
	
	/*!
	 * @brief コンストラクタ ファイルからのデータを予め読む
	 * @param filename 読み込むファイル名を指定する
	*/
	curve(mc_sim::logger& newlogger, std::string filename);
	
	/*!
	 * @brief プロットデータを補間して返す
	 * @param x 独立変数を指定する
	 * @return double 指定された独立変数xに対応する従属変数yの値
	*/
	double itpl_getter(double x);
	
	/*!
	 * @brief 従属変数方データの最大値を返す
	 * @return 従属変数方データの最大値
	*/
	double max();
	
	/*!
	 * @brief データを追加する
	 * @param x 独立変数を指定する
	 * @param y 従属変数を指定する
	*/
	void append(double x, double y);
	
	/*!
	 * @brief データをファイルから読み込む
	 * @param filename ファイル名を指定する
	*/
	void file_input(std::string filename);
	
	/*!
	 * @brief データをファイルに書き出す
	 * @param filename ファイル名を指定する
	*/
	void file_output(std::string filename);
};
