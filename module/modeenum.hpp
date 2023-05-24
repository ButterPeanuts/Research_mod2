/*!
 * @file modeenum.hpp
 * @brief 粒子モードを表す列挙子
 * @author ButterPeanuts
 * @data 2023-05-10
*/
#pragma once

/*!
 * @enum wave_direction
 * @brief 縦波か横波か
*/
enum wave_direction{longitudinal, transverse};

/*!
 * @enum wave_mode
 * @brief 音響モードか光学モードか
*/
enum wave_mode{acoustic, optical};

/*!
 * @enum scatconst
 * @brief 散乱係数を表す
*/
enum scatconst{a, b, chi, xi, c, f};
