#include <iostream>
#include "./src/Simulation2D.h"

using namespace std;

int ngrid = 128;

Simulation2D sim = Simulation2D(ngrid);

void makeScene() {
	double dx = 1.0 / (double)ngrid;
	// 可以给一个初速度好debug
	//sim.sample_drop_rectangle(Vector2d(0.0,-0.58), Vector2d(0.0, 0.0),10000, 0.2);
	sim.sample_drop_rectangle(Vector2d(0.0, 0.3), Vector2d(0.0, -1.0), 10000, 0.2);

	//sim.addSnow(dvec3(0.0, 0.0, 0.0), dvec3(0.0, 0.0, 0.0), 0.2, 1.0, 8000, dvec3(1.0, 1.0, 1.0));
	cout << "makeScene done." << endl;
}

#define ITER 10
#define FRAME_NUM 60

int main() {
	makeScene();
	int frame = 0;
	while (frame < 200) {
		for(int i = 0 ; i <ITER * 20; i ++ )
			sim.update_elastic_pic(1.0 / FRAME_NUM / (ITER * 20));
		sim.render_elastic(frame);
		frame++;
		cout << "frame:" << frame << endl;
	}

	return 0;
}