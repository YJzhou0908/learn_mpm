#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "./src/Simulation2D.h"

using namespace std;

int ngrid = 200;

Simulation2D sim = Simulation2D(ngrid);

void makeScene() {
	double dx = 1.0 / (double)ngrid;
	// ���Ը�һ�����ٶȺ�debug
	sim.sample_drop_snow(Vector2d(0.0, 0.0), Vector2d(0.0, -2.0), 10000, 0.1);

	//sim.addSnow(dvec3(0.0, 0.0, 0.0), dvec3(0.0, 0.0, 0.0), 0.2, 1.0, 8000, dvec3(1.0, 1.0, 1.0));
	cout << "makeScene done." << endl;
}

#define ITER 10
#define FRAME_NUM 60

int main() {
	makeScene();
	int frame = 0;
	while (frame < 200) {
		sim.render_snow(frame);
		for (int i = 0; i < ITER * 10; i++)
			sim.update_snow_apic(1.0 / FRAME_NUM / (ITER * 30));

		frame++;
		cout << "frame:" << frame << endl;
	}

	return 0;
}