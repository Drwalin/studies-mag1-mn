
#include <cstdio>
#include <cmath>

struct Pair {
	float x, y;
};

Pair Solve(Pair p, float dx) {
	Pair r;
	r.x = p.x + dx;
	float dfdx = p.y*p.y + 1;
	r.y = p.y + dfdx*dx;
	return r;
}

void PrintSolutions(float x0, float y0, float dx, float maxx) {
	printf("y(x_0 = %f) = %f\n", x0, y0);
	printf(" dx = %f\n", dx);
	Pair p{x0, y0};
	for(; p.x < maxx+0.00001;) {
		printf(" y( %f ) = %f\n", p.x, p.y);
		p = Solve(p, dx);
	}
	printf("\n\n");
}

int main() {
	float x0 = 0, y0 = 0;
	PrintSolutions(x0, y0, 0.2, 1);
	PrintSolutions(x0, y0, 0.1, 1);
	PrintSolutions(x0, y0, 0.05, 1);
	return 0;
}

