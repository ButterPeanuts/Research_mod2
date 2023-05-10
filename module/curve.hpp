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

/**
 * @brief 二次元グラフ上にプロットできるデータを示すオブジェクトのクラス
*/
class curve{
	public:
	/*! プロットデータを表す */
	std::vector<std::pair<double, double>> table;
	
	/*!
	 * @brief コンストラクタ
	*/
	curve();
	
	/*!
	 * @brief コンストラクタ ファイルからのデータを予め読む
	 * @param filename 読み込むファイル名を指定する
	*/
	curve(std::string filename);
	
	/*!
	 * @brief プロットデータを補間して返す
	 * @param x 独立変数を指定する
	*/
	double itpl_getter(double x);
	
	/*!
	 * @brief 従属変数方データの最大値を返す
	*/
	double max();
	
	/*!
	 * @brief データを追加する
	 * @param x 独立変数を指定する
	 * @param y 従属変数を指定する
	*/
	double append(double x, double y);
	
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
