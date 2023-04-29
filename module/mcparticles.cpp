#include "mcparticles.hpp"
#include "physconst.hpp"
#include "massconst.hpp"
#include<vector>
#include<cmath>
#include<numbers>
#include<iostream>

using namespace mc_particles;


MCParticles::MCParticles(double Energy,double Temperature){
	//�������u(���x����)����
	//https://qiita.com/aa_debdeb/items/e416ae8a018692fc07eb ���Q�Ƃ̂���
	this->velocity_pointing = std::vector<double>(MCParticles::dimension, 0);
	this->Elastic_scattering();
	
	//�������omega,band����
	//��Ԗ��x�̏���MCParticle���Ƃɂǂ̒l��ێ����邩�ς�肻���Ȃ̂�
	//�����o�[��ǉ����Ă����K�v������
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	std::uniform_real_distribution<> randx(1.0e+12, std::max((*(massconst::Si_DOS_LA.end() - 1))[0],(*(massconst::Si_DOS_TA.end() - 1))[0]));
	//���p�@�̂��߂ɍő�l���o���Ă���Ƃ���
	//���ʂ�����
	//���O�Ƀe�[�u���`���Ŏ����Ă����ׂ�
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	double maxdis = 0;
	for (auto i = massconst::Si_DOS_LA.begin(); i < massconst::Si_DOS_LA.end(); i++) {
		double P = (*i)[1] * physconst::dirac * (*i)[0] / Energy / (exp(physconst::dirac * (*i)[0] / physconst::boltzmann / Temperature) - 1);
		if (maxdis > P)maxdis = P;
	}
	
	//�����P���L���ł̓_������
	//�Ȃɂ�DOS_interpolation������̂�
	//�Ă�����DOS_interpolation�̓o���h�N���X�Ɉڊǂł�����
	std::uniform_real_distribution<> randf(0, maxdis);
	std::uniform_int_distribution<> randp(0, 2);
	for (;;) {
		double xr = randx(physconst::mtrand);
		double fr = randf(physconst::mtrand);
		int pr = randp(physconst::mtrand);
		auto P = [=](double omega, int p) {
			if (p == 2) {
				return massconst::DOS_interpolation(massconst::Si_DOS_LA, omega) * omega * physconst::dirac / Energy / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
			}
			else {
				return massconst::DOS_interpolation(massconst::Si_DOS_TA, omega) * omega * physconst::dirac / Energy / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
			}
		};
		if (fr <= P(xr, pr)) {
			this->angular_frequency = xr;
			this->bandnum = pr;
			break;
		}
	}
	
	//�������r����
	this->position[0] = 0;
	this->position[1] = 0;
	this->position[2] = 0;
}

void MCParticles::Nextstep(double dt) {
	//������o���h�N���X�ɓ��������ق����ǂ�����
	//(https://www.notion.so/MCParticle-761aeb843f10432d81f7b255734779fe?pvs=4)
	double velocity = massconst::Si_group_velocity(angular_frequency,bandnum);
	
	for (int i = 0; i < MCParticles::dimension; i++) {
		position[i] += dt * velocity_pointing[i] * velocity;
	}
}

void MCParticles::Boundary_Scatter_B(double max_x, double max_y, double max_z) {
	//�����͋C�ɂȂ�Ƃ��Ɍ���and����t���ł����̂�
	if ((this->position)[1] < 0 || max_y < this->position[1]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[1] * this->velocity_pointing[1]);
		velocity_pointing[0] /= sin_oldtheta;
		velocity_pointing[2] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[1] = (this->position[1] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
	}
	if ((this->position)[0] < 0 || max_x < this->position[0]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[0] * this->velocity_pointing[0]);
		velocity_pointing[1] /= sin_oldtheta;
		velocity_pointing[2] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[0] = (this->position[0] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[1] *= sin_newtheta;
		velocity_pointing[2] *= sin_newtheta;
	}
	if ((this->position)[2] < 0 || max_z < this->position[2]) {
		std::uniform_real_distribution<> randR(0, 1);
		double sin_oldtheta = std::sqrt(1 - this->velocity_pointing[2] * this->velocity_pointing[2]);
		velocity_pointing[0] /= sin_oldtheta;
		velocity_pointing[1] /= sin_oldtheta;
		double cos_newtheta = std::sqrt(randR(physconst::mtrand));
		velocity_pointing[2] = (this->position[2] < 0 ? 1 : -1) * cos_newtheta;
		double sin_newtheta = std::sqrt(1 - cos_newtheta * cos_newtheta);
		velocity_pointing[0] *= sin_newtheta;
		velocity_pointing[1] *= sin_newtheta;
	}
}

void MCParticles::Scatter(double Temperature,double dt,double min_structure) {
    //���E�U��B���K�v
    //�t�H�m���t�H�m���U��
    //�o���h�ԍ��Ȃ�Ƃ����Ȃ��Ƃ�
	std::uniform_real_distribution<> randx(0, 1);
    std::uniform_real_distribution<> randcosth(-1, 1);

    //�t�H�m������
    double Pu;
    if (this->bandnum == 0 || this->bandnum == 1) {
        Pu = massconst::Si_scatter_ATA * pow(this->angular_frequency,massconst::Si_scatter_chiTA) * pow(Temperature,massconst::Si_scatter_xiTA) * std::exp(massconst::Si_scatter_BTA / (-Temperature));
    }
    else {
        Pu = massconst::Si_scatter_ALA * pow(this->angular_frequency,massconst::Si_scatter_chiLA) * pow(Temperature,massconst::Si_scatter_xiLA) * std::exp(massconst::Si_scatter_BLA / (-Temperature));
    }
    if (randx(physconst::mtrand) <= (1 - exp(-dt * Pu))) {
		std::uniform_real_distribution<> randx(std::min((*(massconst::Si_DOS_LA.begin() + 1))[0],(*(massconst::Si_DOS_TA.begin() + 1))[0]), std::max((*(massconst::Si_DOS_LA.end() - 1))[0],(*(massconst::Si_DOS_TA.end() - 1))[0]));
		double maxdis = 0;
		for (auto i = massconst::Si_DOS_LA.begin(); i < massconst::Si_DOS_LA.end(); i++) {
			double P = Pu * (*i)[1] * physconst::dirac * (*i)[0] / (physconst::dirac * this->angular_frequency)/ (exp(physconst::dirac * (*i)[0] / physconst::boltzmann / Temperature) - 1);
			if (maxdis > P)maxdis = P;
		}

		std::uniform_real_distribution<> randf(0, maxdis);
		std::uniform_int_distribution<> randp(0, 2);
		for (;;) {
			double xr = randx(physconst::mtrand);
			double fr = randf(physconst::mtrand);
			int pr = randp(physconst::mtrand);
			auto P = [&](double omega, int p) {
				if (p == 2) {
					return Pu * massconst::DOS_interpolation(massconst::Si_DOS_LA, omega) * omega * physconst::dirac / (physconst::dirac * this->angular_frequency) / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
				}
				else {
					return Pu * massconst::DOS_interpolation(massconst::Si_DOS_TA, omega) * omega * physconst::dirac / (physconst::dirac * this->angular_frequency) / (exp(physconst::dirac * omega / physconst::boltzmann / Temperature) - 1);
				}
			};
			if (fr <= P(xr, pr)) {
				this->angular_frequency = xr;
				this->bandnum = pr;
				break;
			}
		}
        return;
    }

    //�t�H�m������
    double Pd = massconst::Si_scatter_C * pow(angular_frequency, 4);
    if (randx(physconst::mtrand) <= (1 - exp(-dt * Pd))) {
        this->Elastic_scattering();
        return;
    }

    //�t�H�m�����EA
    //F,L�͌��
    double Pb = massconst::Si_group_velocity(angular_frequency, bandnum) * min_structure * 0.55;
    if (randx(physconst::mtrand) <= (1 - exp(-dt * Pb))) {
        this->Elastic_scattering();
        return;
    }
}

//�e���U��(�����ω��Ȃ�,���x�x�N�g�������ω�)
void MCParticles::Elastic_scattering() {
	std::uniform_real_distribution<> randcosth(-1, 1);
	double costh = randcosth(physconst::mtrand);
	double phi = randcosth(physconst::mtrand) * std::numbers::pi;
	double sinth = std::sqrt(1 - costh * costh);
	this->velocity_pointing[0] = sinth * std::cos(phi);
	this->velocity_pointing[1] = sinth * std::sin(phi);
	this->velocity_pointing[2] = costh;
	return;
}
