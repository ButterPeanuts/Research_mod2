#include<random>
#include "modeenum.hpp"
#include "band_obj.hpp"
#include <logger.hpp>
#include <scatconst.hpp>

band_obj::band_obj(std::shared_ptr<mc_sim::logger>& newlogger, curve& dos, curve& gvelocity, curve& domcpmax, wave_direction direction, wave_mode mode, scatconst& bands_scatconst) : dos(dos), gvelocity(gvelocity), domcpmax(domcpmax), bands_scatconst(bands_scatconst), logger(newlogger){
	this->directions = direction;
	this->mode = mode;
	this->dos_max = this->dos.max();
}

double band_obj::dos_getter(double omega){
	return this->dos.itpl_getter(omega);
}

std::uniform_real_distribution<> band_obj::dos_omega_distribution_getter(){
	return std::uniform_real_distribution<>(this->dos.left_edge(), this->dos.right_edge());
}

std::uniform_real_distribution<> band_obj::dos_distribution_getter(){
	return std::uniform_real_distribution<>(0, this->dos_max);
}

std::uniform_real_distribution<> band_obj::domcp_distribution_getter(double t){
	return std::uniform_real_distribution<>(0, this->domcpmax.itpl_getter(t));
}

double band_obj::gvelocity_getter(double omega){
	return this->gvelocity.itpl_getter(omega);
}

double band_obj::dos_leftedge(){
	return this->dos.left_edge();
}

double band_obj::dos_rightedge(){
	return this->dos.right_edge();
}

/* double band_obj::domcp_getter(double omega){ */
/* 	return this->domcp.itpl_getter(omega); */
/* } */

/* std::uniform_real_distribution<> band_obj::domcp_omega_distribution_getter(){ */
/* 	return std::uniform_real_distribution<>(this->domcp.left_edge(), this->domcp.right_edge()); */
/* } */

/* std::uniform_real_distribution<> band_obj::domcp_distribution_getter(){ */
/* 	return std::uniform_real_distribution<>(0, this->domcp_max); */
/* } */

wave_direction band_obj::directions_getter(){
	return this->directions;
}

wave_mode band_obj::mode_getter(){
	return this->mode;
}

double band_obj::a(){
	return this->bands_scatconst.a();
}

double band_obj::b(){
	return this->bands_scatconst.b();
}

double band_obj::chi(){
	return this->bands_scatconst.chi();
}

double band_obj::xi(){
	return this->bands_scatconst.xi();
}

double band_obj::c(){
	return this->bands_scatconst.c();
}

double band_obj::f(){
	return this->bands_scatconst.f();
}

double band_obj::scatprob_u(double omega, double t, double dt){
	return -std::expm1(-this->bands_scatconst.tau_u_inv(omega, t) * dt);
}

double band_obj::scatprob_d(double omega, double dt){
	return -std::expm1(-this->bands_scatconst.tau_d_inv(omega) * dt);
}

double band_obj::scatprob_b(double gvelocity, double l, double dt){
	return -std::expm1(-this->bands_scatconst.tau_b_inv(gvelocity, l) * dt);
}

double band_obj::mintau(double tmax, double l){
	auto& scat = this->bands_scatconst;
	double omegamax = this->dos_rightedge();
	double gvelocitymax = this->gvelocity.max();
	return 1 / (scat.tau_u_inv(omegamax, tmax) + scat.tau_d_inv(omegamax) + scat.tau_b_inv(gvelocitymax, l));
}
