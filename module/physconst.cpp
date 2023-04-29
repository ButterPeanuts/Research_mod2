#include"physconst.hpp"
#include<random>
static std::random_device thrand;

std::mt19937_64 physconst::mtrand(thrand());
double physconst::vonNeumann_rejection(double (*f)(double),double xi, double xs, double fm){
	std::uniform_real_distribution<> randx(xi, xs);
	std::uniform_real_distribution<> randf(0, fm);
	for (;;) {
		double xr = randx(physconst::mtrand);
		double fr = randf(physconst::mtrand);
		if (fr <= f(xr))return xr;
	}
}
