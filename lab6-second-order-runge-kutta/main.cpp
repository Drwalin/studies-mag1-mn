
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>

using Real = double;
struct Pair {
	Real x, y;
};
struct Tri {
	Real x, y, yb;
	operator Pair() {
		return {x, y};
	}
};
Real fabs(Real v) {
	if(v < 0) return -v; else return v;
}


Real f(Tri p) {
	return -sin(p.y) -0.5*p.yb + 0.9 * cos(0.45*p.x);
}
Real f(Real x, Real y, Real yb) {
	return f({x,y,yb});
}

Tri SimpleEuler(Tri arg, Real dx) {
	Tri r = arg;
	r.x += dx;
	r.y += r.yb*dx;
	r.yb += f(arg)*dx;
	return r;
}

Tri MoreImprovedEulerByMarekZalewski(Tri arg, Real dx) {
	Tri r = arg;
	r.x += dx;
	r.yb += (f(arg) + f(SimpleEuler(arg, dx)))*dx/2;
	r.y += (arg.yb + r.yb)*dx/2;
	return r;
}

Tri ImprovedEuler(Tri arg, Real dx) {
	Tri r = arg;
	r.x += dx;
	r.y += (arg.yb*2 + f(arg)*dx)*dx/2;
	r.yb += (f(arg) + f(SimpleEuler(arg, dx)))*dx/2;
	return r;
}

// Tri RungeKuttaReferenceImplementation(Tri arg, Real dx) {
// 	Real f0, f1, f2, f3;
// 	Tri r = arg;
// 	r.x += dx;
// 	
// 	Real f00 = arg.yb;
// 	Real f01 = f(arg);
// 	Real f10 = arg.yb + f01*dx/2;
// 	Real f11 = f(arg.x+dx/2, arg.y+f00*dx/2, arg.yb+f01*dx/2);
// 	Real f20 = arg.yb + f11*dx/2;
// 	Real f21 = f(arg.x+dx/2, arg.y+f10*dx/2, arg.yb+f11*dx/2);
// 	Real f30 = arg.yb+f21*dx;
// 	Real f31 = f(arg.x+dx, arg.y+f20*dx, arg.yb+f21*dx);
// 	
// 	r.y += (f00 + 2*f10 + 2*f20 + f30)*dx/6;
// 	r.yb += (f01 + 2*f11 + 2*f21 + f31)*dx/6;
// 	return r;
// }

Tri RungeKutta(Tri arg, Real dx) {
	Tri r = arg;
	r.x += dx;
	
	Real f00 = arg.yb;
	Real f01 = f(arg);
	
	Real f10 = arg.yb + f01*dx/2;
	Real f11 = f(arg.x+dx/2, arg.y+f00*dx/2, arg.yb+f01*dx/2);
	
	Real f20 = arg.yb + f11*dx/2;
	Real f21 = f(arg.x+dx/2, arg.y+f10*dx/2, arg.yb+f11*dx/2);
	
	Real f30 = arg.yb + f21*dx;
	Real f31 = f(arg.x+dx, arg.y+f20*dx, arg.yb+f21*dx);
	
	r.yb += (f01 + 2*f11 + 2*f21 + f31)*dx/6;
	r.y  += (f00 + 2*f10 + 2*f20 + f30)*dx/6;
	return r;
}



std::vector<Tri> SolveSecondOrder(Tri p0, Real dx, Real maxx,
		Tri(*func)(Tri,Real)) {
	std::vector<Tri> ret;
	Tri p = p0;
	ret.emplace_back(p0);
	do {
		p = func(p, dx);
		ret.emplace_back(p);
	} while(p.x+0.01 <= maxx);
	return ret;
}

using FFtypeF = Tri(*)(Tri,Real);
std::vector<std::vector<Tri>> SolveSecondOrderMethods(Tri p0, Real dx, Real maxx,
		std::vector<FFtypeF> func) {
	std::vector<std::vector<Tri>> ret;
	for(int i=0; i<func.size(); ++i) {
		ret.emplace_back(
				SolveSecondOrder(p0, dx, maxx, func[i])
			);
	}
	return ret;
}

void PrintSolutions(const std::vector<std::vector<Tri>>& sol) {
	for(int i=0; i<sol.size(); ++i) {
		for(const auto& v : sol[i]) {
			printf(" f(%f) = \t%f\n", v.x, v.y);
		}
		printf("\n\n");
	}
}

int main() {
	auto sol = SolveSecondOrderMethods({0,1.4,0}, 0.1, 3.0, {
// 		RungeKuttaReferenceImplementation,
		RungeKutta});
	
	PrintSolutions(sol);
	
	return 0;
}

