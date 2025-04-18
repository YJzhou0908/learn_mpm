#include "common.h"
#include "Frame.h"
#include "Particle2D.h"
#include "Grid2D.h"

using namespace std;
using vec = Vector2d;
class Simulation2D {
public:
	uint32_t grid_num = 128;
	double dx = 1.0 / grid_num;
	Grid2D mGrid;
	vector<Particle2D> mParticles;
	Frame mFrame;

	const double mE = 1.4e5;// young's Modulus // fem的stiffness matrix的condition number， CLF condition，仿真 1e4 0.5 ,仿真部署1/60/5
	const double mNu = 0.2; // Possion ratio
	const double hardening = 2.0; // 2.0 和 7.5e-2
	const double flipAlpha = 0.95;
	const double ground = 0.0;
	const double mCompression = 5.5e-2;
	const double mStrech = 7.5e-3;
	const Vector2d mGravity = Vector2d(0.0, -9.8);

	const double mu0 = mE / (2 * (1 + mNu));
	const double lambda0 = (mE * mNu) / ((1 + mNu) * (1 - 2 * mNu));


	Simulation2D(uint32_t ngrid) {
		dx = 1.0 / ngrid;
		mGrid = Grid2D(ngrid);
		grid_num = ngrid;
		mFrame = Frame(600, 600, 3);
	}

	void sample_rotate_rectangle(Vector2d&c, int n, double hight, double weight, double omega) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dist(-hight / 2, hight / 2);
		int i = 0;
		double p_mass = weight * hight * hight * 0.93 / n;
		while (i < n) {
			double x = dist(gen);
			double y = dist(gen);
			Vector2d vel = omega * Vector2d(-y, x);
			mParticles.push_back(Particle2D(Vector2d(c(0) + x, c(1) + y), vel, p_mass, 0.0, 0));
			if (c(1) + y > 1 || c(1) + y < -1)
				cout << y << endl;
			i++;
		}
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
		}
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					particle.mDensity += nodeMass * weight;
				}
			}
			particle.mDensity /= pow(mGrid.mCellSize, 3);
			particle.mVol = particle.mMass / particle.mDensity;
		}

		mGrid.clear();
	
	}

	void sample_drop_snow(Vector2d& c, Vector2d& vel, int n, double radius, double p_mass = 0.000000248, double p_vol = 0.0){
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dist(-1.0, 1.0);
			std::uniform_real_distribution<> radiusDist(0.0, 1.0);
			
			p_mass = M_PI * radius * radius * 0.4 / n;
			int i = 0;
			while (i < n) {
				double x = dist(gen);
				double y = dist(gen);
				double lenSq = x * x + y * y;
				if (lenSq <= 1.0 && lenSq > 0.0) {
					double len = std::sqrt(lenSq);
					x /= len;
					y /= len;
					double scale = std::cbrt(radiusDist(gen)) * radius;
					mParticles.push_back(Particle2D(Vector2d(c(0) + x * scale, c(1) + y * scale), vel, p_mass, p_vol, 0));
					i++;
				}
			}
			
			for (int i = 0; i < mParticles.size(); i++) {
				vec& pos = mParticles[i].mPosition;
				auto& particle = mParticles[i];
				int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
				int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
				vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {
						int grid_x = c_x + r;
						int grid_y = c_y + c;
						// delta x , 注意：单位一定是网格！！
						vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
						double weight = computeWeightQuad(dx); // 
						vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
						mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
					}
				}
			}
			for (int i = 0; i < mParticles.size(); i++) {
				vec& pos = mParticles[i].mPosition;
				auto& particle = mParticles[i];
				vec& velocity = particle.mVelocity;
				int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
				int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
				vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
				//pMomentum += mParticles[i].mMass * velocity.norm();
				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {
						int grid_x = c_x + r;
						int grid_y = c_y + c;
						double nodeMass = mGrid.getMass(grid_x, grid_y);
						vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
						double weight = computeWeightQuad(dx);
						particle.mDensity += nodeMass * weight;
					}
				}
				particle.mDensity /= pow(mGrid.mCellSize, 3);
				particle.mVol = particle.mMass / particle.mDensity;
			}

			mGrid.clear();
		
	
	}

	void sample_drop_rectangle(Vector2d& c, Vector2d& vel, int n, double edge, double p_mass = 0.000000248, double p_vol = 0.0) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dist(-edge / 2, edge / 2);
		int i = 0;
		p_mass = edge * edge * 0.93 / n;
		while (i < n) {
			double x = dist(gen);
			double y = dist(gen);
			mParticles.push_back(Particle2D(Vector2d(c(0) + x, c(1) + y), vel, p_mass, p_vol, 0));
			i++;
		}
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
		}
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					particle.mDensity += nodeMass * weight;
				}
			}
			particle.mDensity /= pow(mGrid.mCellSize, 3);
			particle.mVol = particle.mMass / particle.mDensity;
		}

		mGrid.clear();
	}
	
	void sample_drop_rectangle_material(Vector2d& c, Vector2d& vel, int n, double edge, uint32_t material_type) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dist(-edge / 2, edge / 2);
		
		double p_mass =  edge* edge / n;
		if (material_type == SNOW) {
			p_mass *= 0.93;
		}
		else if (material_type == ELASITC) {
			p_mass *= 0.4;
		}
		int i = 0;
		int pBegin = mParticles.size();
		int pEnd = pBegin + n;
		while (i < n) {
			double x = dist(gen);
			double y = dist(gen);
			mParticles.push_back(Particle2D(Vector2d(c(0) + x, c(1) + y), vel, p_mass, 0.0, 0, material_type));
			i++;
		}

		for (int i = pBegin; i < pEnd; i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
		}
		for (int i = pBegin; i < pEnd; i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					particle.mDensity += nodeMass * weight;
				}
			}
			particle.mDensity /= pow(mGrid.mCellSize, 3);
			particle.mVol = particle.mMass / particle.mDensity;
		}

		mGrid.clear();
	}

	void sample_drop_ball_material(Vector2d& c, Vector2d& vel, int n, double radius, uint32_t material_type) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dist(-1.0, 1.0);
		std::uniform_real_distribution<> radiusDist(0.0, 1.0);

		double p_mass = M_PI * radius * radius  / n;
		if (material_type == SNOW) {
			p_mass *= 0.93;
		}
		else if (material_type == ELASITC) {
			p_mass *= 0.4;
		}

		int pBegin = mParticles.size();
		int pEnd = pBegin + n;

		int i = 0;
		while (i < n) {
			double x = dist(gen);
			double y = dist(gen);
			double lenSq = x * x + y * y;
			if (lenSq <= 1.0 && lenSq > 0.0) {
				double len = std::sqrt(lenSq);
				x /= len;
				y /= len;
				double scale = std::cbrt(radiusDist(gen)) * radius;
				mParticles.push_back(Particle2D(Vector2d(c(0) + x * scale, c(1) + y * scale), vel, p_mass, 0.0, 0, material_type));
				i++;
			}
		}

		for (int i = pBegin; i < pEnd; i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
		}
		for (int i = pBegin; i < pEnd; i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					particle.mDensity += nodeMass * weight;
				}
			}
			particle.mDensity /= pow(mGrid.mCellSize, 3);
			particle.mVol = particle.mMass / particle.mDensity;
		}

		mGrid.clear();

	}

	void update_rigid_apic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

		double pMomentum = 0.0;
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, (particle.mMass * velocity + particle.mMass * particle.mCp * dx * mGrid.mCellSize) * weight / nodeMass);

				}
			}
		}
		
		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			particle.mCp = Matrix2d::Zero();
			vec V_apic = vec(0.0, 0.0);
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_apic += weight * node.mVelocity;
					particle.mCp += 4 * weight / mGrid.mCellSize * node.mVelocity * dx.transpose() ;

				}
			}
			particle.mVelocity = V_apic;
		}
		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}
		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	}

	void update_rigid_pic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

		double pMomentum = 0.0;
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, particle.mMass * velocity * weight / nodeMass);

				}
			}
		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_pic = vec(0.0, 0.0);
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_pic += weight * node.mVelocity;

				}
			}
			particle.mVelocity = V_pic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	
	}

	void update_elastic_pic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

		double pMomentum = 0.0;
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); 
					mGrid.addVelocity(grid_x, grid_y, particle.mMass * velocity  * weight / nodeMass);
				}
			}
		}
		/*
		double gMomentum = 0.0;
		for (int i = 0; i < mGrid.mNode.size(); i++) {
			gMomentum += mGrid.mNode[i].mMass * mGrid.mNode[i].mVelocity.norm();
		}
		*/
		
		for (int i = 0; i < mParticles.size(); i++)
		{
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			auto& Fe = particle.mFe;
			double Je = particle.mDetFe;
			auto FeT = Fe.transpose();
			double mu = mu0;
			double lambda = lambda0;
			JacobiSVD<Eigen::Matrix2d> svd(Fe,
				Eigen::ComputeFullU |
				Eigen::ComputeFullV);
			const Matrix2d& U = svd.matrixU();
			const Matrix2d& V = svd.matrixV();
			const Matrix2d& S = svd.singularValues().asDiagonal();
			const Matrix2d Re = U * V.transpose(); // 极分解 
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					//double weight = computeWeight(dx); // 
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					// f = -Vp / J * (2*mu*(F-R)*FT + lambda * (J - 1) * J  * I) * dw;
					vec force = //  -particle.mVol * lambda * (Je - 1.0) * Je * dmat3(1.0) * dw;
						-particle.mVol * (2.0 * mu * (Fe - Re) * FeT + lambda * (Je - 1.0) * Je * Matrix2d::Identity()) * dw;
					mGrid.addForce(grid_x, grid_y, force);
				}
			}

		}

		for (int i = 0; i < mGrid.mNode.size(); i++) {
			mGrid.mNode[i].mForce += mGrid.mNode[i].mMass * mGravity;
			if (mGrid.mNode[i].mMass == 0.0) {
				mGrid.mNode[i].mVelocity(0) = 0.0;
				mGrid.mNode[i].mVelocity(1) = 0.0;
			}
			else {
				mGrid.mNode[i].mVelocity += mGrid.mNode[i].mForce / mGrid.mNode[i].mMass * dt;
			}
		}

		for (int i = 0; i < grid_num + 1; i++) {
			mGrid.mNode[i].mVelocity(1) = 0.0;
		}
		
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			Matrix2d gradV = Matrix2d::Zero();
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					vec vel = mGrid.getVelocity(grid_x, grid_y);
					gradV += vel * dw.transpose();
				}
			}
			particle.mFe = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
			particle.mDetFe = particle.mFe.determinant();
			//if (1.0 - particle.mDetFe > 1e-5)
			//	cout << particle.mDetFe << endl;
		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_pic = vec(0.0, 0.0);
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_pic += weight * node.mVelocity;
					
				}
			}
			particle.mVelocity = V_pic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	}

	void update_elastic_apic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

		double pMomentum = 0.0;
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, (particle.mMass * velocity + particle.mMass * particle.mCp * dx * mGrid.mCellSize)* weight / nodeMass);

				}
			}
		}

		for (int i = 0; i < mParticles.size(); i++)
		{
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			auto& Fe = particle.mFe;
			double Je = particle.mDetFe;
			auto FeT = Fe.transpose();
			double mu = mu0;
			double lambda = lambda0;
			JacobiSVD<Eigen::Matrix2d> svd(Fe,
				Eigen::ComputeFullU |
				Eigen::ComputeFullV);
			const Matrix2d& U = svd.matrixU();
			const Matrix2d& V = svd.matrixV();
			const Matrix2d& S = svd.singularValues().asDiagonal();
			const Matrix2d Re = U * V.transpose(); // 极分解 
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					// f = -Vp / J * (2*mu*(F-R)*FT + lambda * (J - 1) * J  * I) * dw;
					vec force = //  -particle.mVol * lambda * (Je - 1.0) * Je * dmat3(1.0) * dw;
						-particle.mVol * (2.0 * mu * (Fe - Re) * FeT + lambda * (Je - 1.0) * Je * Matrix2d::Identity()) * dw;
					mGrid.addForce(grid_x, grid_y, force);
				}
			}

		}

		for (int i = 0; i < mGrid.mNode.size(); i++) {
			//mGrid.mNode[i].mForce += mGrid.mNode[i].mMass * mGravity;
			if (mGrid.mNode[i].mMass == 0.0) {
				mGrid.mNode[i].mVelocity(0) = 0.0;
				mGrid.mNode[i].mVelocity(1) = 0.0;
			}
			else {
				mGrid.mNode[i].mVelocity += mGrid.mNode[i].mForce / mGrid.mNode[i].mMass * dt;
			}
		}

		for (int i = 0; i < grid_num + 1; i++) {
			mGrid.mNode[i].mVelocity(1) = 0.0;
		}

		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			Matrix2d gradV = Matrix2d::Zero();
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					vec vel = mGrid.getVelocity(grid_x, grid_y);
					gradV += vel * dw.transpose();
				}
			}
			particle.mFe = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
			particle.mDetFe = particle.mFe.determinant();
		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_apic = vec(0.0, 0.0);
			particle.mCp = Matrix2d::Zero();
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_apic += weight * node.mVelocity;
					particle.mCp += 4 * weight / mGrid.mCellSize * node.mVelocity * dx.transpose();
				}
			}
			particle.mVelocity = V_apic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	}

	void update_snow_apic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

		double pMomentum = 0.0;
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, (particle.mMass * velocity + particle.mMass * particle.mCp * dx * mGrid.mCellSize) * weight / nodeMass);

				}
			}
		}

		for (int i = 0; i < mParticles.size(); i++)
		{
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			auto& Fe = particle.mFe;
			auto& Fp = particle.mFp;
			double Jp = particle.mDetFp;
			double Je = particle.mDetFe;
			double J = Jp * Je;
			auto FeT = Fe.transpose();
			double mu = mu0 * std::exp(hardening * (1.0 - Jp));
			double lambda = lambda0 * std::exp(hardening * (1.0 - Jp));
			JacobiSVD<Eigen::Matrix2d> svd(Fe,
				Eigen::ComputeFullU |
				Eigen::ComputeFullV);
			const Matrix2d& U = svd.matrixU();
			const Matrix2d& V = svd.matrixV();
			const Matrix2d& S = svd.singularValues().asDiagonal();
			const Matrix2d Re = U * V.transpose(); // 极分解 
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					// f = -Vp / J * (2*mu*(F-R)*FT + lambda * (J - 1) * J  * I) * dw;
					vec force = //  -particle.mVol * lambda * (Je - 1.0) * Je * dmat3(1.0) * dw;
						-particle.mVol * (2.0 * mu * (Fe - Re) * FeT + lambda * (Je - 1.0) * Je * Matrix2d::Identity()) * dw;
					mGrid.addForce(grid_x, grid_y, force);
				}
			}

		}

		for (int i = 0; i < mGrid.mNode.size(); i++) {
			mGrid.mNode[i].mForce += mGrid.mNode[i].mMass * mGravity;
			mGrid.mNode[i].mOldVelocity = mGrid.mNode[i].mVelocity;
			if (mGrid.mNode[i].mMass == 0.0) {
				mGrid.mNode[i].mVelocity(0) = 0.0;
				mGrid.mNode[i].mVelocity(1) = 0.0;
			}
			else {
				mGrid.mNode[i].mVelocity += mGrid.mNode[i].mForce / mGrid.mNode[i].mMass * dt;
			}
		}

		for (int i = 0; i < grid_num + 1; i++) {
			mGrid.mNode[i].mVelocity(1) = 0.0;
		}

		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			Matrix2d gradV = Matrix2d::Zero();
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					vec vel = mGrid.getVelocity(grid_x, grid_y);
					gradV += vel * dw.transpose();
				}
			}
			Matrix2d Fe_head = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
			Matrix2d F_next = Fe_head * particle.mFp;

			// clamp Fe_head
			Eigen::JacobiSVD<Eigen::Matrix2d> svd(
				Fe_head, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Vector2d sigma = svd.singularValues();
			double a = 1 - mCompression;
			double b = 1 + mStrech;
			sigma(0) = std::clamp(sigma(0), a, b);
			sigma(1) = std::clamp(sigma(1), a, b);
			particle.mFe = svd.matrixU() * sigma.asDiagonal() * svd.matrixV().transpose();
			particle.mFp = svd.matrixV() * sigma.cwiseInverse().asDiagonal() * svd.matrixU().transpose() * F_next;
			particle.mDetFe = particle.mFe.determinant();
			particle.mDetFp = particle.mFp.determinant();

		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_apic = vec(0.0, 0.0);
			particle.mCp = Matrix2d::Zero();
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_apic += weight * node.mVelocity;
					particle.mCp += 4 * weight / mGrid.mCellSize * node.mVelocity * dx.transpose();
				}
			}
			particle.mVelocity = V_apic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();


	}

	void render_rigid_pic(int frameIdx) {
		std::vector<Vector2d> pixels;
		std::vector<Vector3d> color;
		for (int i = 0; i < mParticles.size(); i++) {
			Vector2d ndcPos = mParticles[i].mPosition; //position 一定是 -1,1 之间的数
			pixels.push_back(ndcPos);
		}

		mFrame.clear();
		mFrame.resize();
		// render ndc position
		for (int i = 0; i < pixels.size(); i++) {
			int pixel_x = static_cast<int>((pixels[i](0) + 1.0) * 0.5 * (mFrame.mWidth));
			int pixel_y = static_cast<int>(((pixels[i](1) + 1.0) * 0.5) * (mFrame.mHeight));
			for (int r = -1; r < 2; r++) {
				for (int c = -1; c < 2; c++) {
					mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 1);
				}
			}
		}
		char buff[100];
		snprintf(buff, sizeof(buff), "rigid_pic%03d.jpg", frameIdx);
		std::string filename = buff;
		mFrame.outputImg(filename);

	}

	void render_rigid_apic(int frameIdx) {
		std::vector<Vector2d> pixels;
		std::vector<Vector3d> color;
		for (int i = 0; i < mParticles.size(); i++) {
			Vector2d ndcPos = mParticles[i].mPosition; //position 一定是 -1,1 之间的数
			pixels.push_back(ndcPos);
		}

		mFrame.clear();
		mFrame.resize();
		// render ndc position
		for (int i = 0; i < pixels.size(); i++) {
			int pixel_x = static_cast<int>((pixels[i](0) + 1.0) * 0.5 * (mFrame.mWidth));
			int pixel_y = static_cast<int>(((pixels[i](1) + 1.0) * 0.5) * (mFrame.mHeight));
			for (int r = -1; r < 2; r++) {
				for (int c = -1; c < 2; c++) {
					mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 1);
				}
			}
		}
		char buff[100];
		snprintf(buff, sizeof(buff), "rigid_apic%03d.jpg", frameIdx);
		std::string filename = buff;
		mFrame.outputImg(filename);

	}

	void render_elastic (int frameIdx) {
		std::vector<Vector2d> pixels;
		std::vector<Vector3d> color;
		for (int i = 0; i < mParticles.size(); i++) {
			Vector2d ndcPos = mParticles[i].mPosition; //position 一定是 -1,1 之间的数
			pixels.push_back(ndcPos);
		}

		mFrame.clear();
		mFrame.resize();
		// render ndc position
		for (int i = 0; i < pixels.size(); i++) {
			int pixel_x = static_cast<int>((pixels[i](0) + 1.0) * 0.5 * (mFrame.mWidth));
			int pixel_y = static_cast<int>(((pixels[i](1) + 1.0) * 0.5) * (mFrame.mHeight));
			for (int r = -1; r < 2; r++) {
				for (int c = -1; c < 2; c++) {
					mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 1);
				}
			}
		}
		char buff[100];
		snprintf(buff, sizeof(buff), "elastic%03d.jpg", frameIdx);
		std::string filename = buff;
		mFrame.outputImg(filename);

	}

	void render_snow(int frameIdx) {
		std::vector<Vector2d> pixels;
		std::vector<Vector3d> color;
		for (int i = 0; i < mParticles.size(); i++) {
			Vector2d ndcPos = mParticles[i].mPosition; //position 一定是 -1,1 之间的数
			pixels.push_back(ndcPos);
		}

		mFrame.clear();
		mFrame.resize();
		// render ndc position
		for (int i = 0; i < pixels.size(); i++) {
			int pixel_x = static_cast<int>((pixels[i](0) + 1.0) * 0.5 * (mFrame.mWidth));
			int pixel_y = static_cast<int>(((pixels[i](1) + 1.0) * 0.5) * (mFrame.mHeight));
			for (int r = -1; r < 2; r++) {
				for (int c = -1; c < 2; c++) {
					mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 1);
					mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 1);
				}
			}
		}
		char buff[100];
		snprintf(buff, sizeof(buff), "snow%03d.jpg", frameIdx);
		std::string filename = buff;
		mFrame.outputImg(filename);
	}

	void update_interaction_pic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}


		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, (particle.mMass * velocity) * weight / nodeMass);

				}
			}
		}


		for (int i = 0; i < mParticles.size(); i++)
		{
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			auto& Fe = particle.mFe;
			auto& Fp = particle.mFp;
			double Jp = particle.mDetFp;
			double Je = particle.mDetFe;
			double J = Jp * Je;
			auto FeT = Fe.transpose();
			double mu = 0.0;
			double lambda = 0.0;
			if (particle.mMaterials.mType == SNOW) {
				mu = particle.mMaterials.mu0 * std::exp(hardening * (1 - Jp));
				lambda = particle.mMaterials.lambda0 * std::exp(hardening * (1 - Jp));
			}
			else if (particle.mMaterials.mType == ELASITC) {
				mu = particle.mMaterials.mu0;
				lambda = particle.mMaterials.lambda0;
			}

			JacobiSVD<Eigen::Matrix2d> svd(Fe,
				Eigen::ComputeFullU |
				Eigen::ComputeFullV);
			const Matrix2d& U = svd.matrixU();
			const Matrix2d& V = svd.matrixV();
			const Matrix2d& S = svd.singularValues().asDiagonal();
			const Matrix2d Re = U * V.transpose(); // 极分解 
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					// f = -Vp / J * (2*mu*(F-R)*FT + lambda * (J - 1) * J  * I) * dw;
					vec force = //  -particle.mVol * lambda * (Je - 1.0) * Je * dmat3(1.0) * dw;
						-particle.mVol * (2.0 * mu * (Fe - Re) * FeT + lambda * (Je - 1.0) * Je * Matrix2d::Identity()) * dw;
					mGrid.addForce(grid_x, grid_y, force);
				}
			}

		}


		for (int i = 0; i < mGrid.mNode.size(); i++) {
			mGrid.mNode[i].mForce += mGrid.mNode[i].mMass * mGravity;
			mGrid.mNode[i].mOldVelocity = mGrid.mNode[i].mVelocity;
			if (mGrid.mNode[i].mMass == 0.0) {
				mGrid.mNode[i].mVelocity(0) = 0.0;
				mGrid.mNode[i].mVelocity(1) = 0.0;
			}
			else {
				mGrid.mNode[i].mVelocity += mGrid.mNode[i].mForce / mGrid.mNode[i].mMass * dt;
			}
		}

		for (int i = 0; i < grid_num + 1; i++) {
			mGrid.mNode[i].mVelocity(1) = 0.0;
		}

		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			Matrix2d gradV = Matrix2d::Zero();
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					vec vel = mGrid.getVelocity(grid_x, grid_y);
					gradV += vel * dw.transpose();
				}
			}
			if (particle.mMaterials.mType == SNOW) {
				Matrix2d Fe_head = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
				Matrix2d F_next = Fe_head * particle.mFp;

				// clamp Fe_head
				Eigen::JacobiSVD<Eigen::Matrix2d> svd(
					Fe_head, Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::Vector2d sigma = svd.singularValues();
				double a = 1 - mCompression;
				double b = 1 + mStrech;
				sigma(0) = std::clamp(sigma(0), a, b);
				sigma(1) = std::clamp(sigma(1), a, b);
				particle.mFe = svd.matrixU() * sigma.asDiagonal() * svd.matrixV().transpose();
				particle.mFp = svd.matrixV() * sigma.cwiseInverse().asDiagonal() * svd.matrixU().transpose() * F_next;
				particle.mDetFe = particle.mFe.determinant();
				particle.mDetFp = particle.mFp.determinant();
			}
			else if (particle.mMaterials.mType == ELASITC) {
				particle.mFe = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
				particle.mDetFe = particle.mFe.determinant();

			}

		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_pic = vec(0.0, 0.0);
			
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_pic += weight * node.mVelocity;
				
				}
			}
			particle.mVelocity = V_pic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	}

	void update_interaction_apic(double dt) {
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					// delta x , 注意：单位一定是网格！！
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					vec dweight = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					mGrid.addMass(grid_x, grid_y, mParticles[i].mMass * weight);
				}
			}
			//assert(abs(sumWeight) - 1 < 1e-3); // 检查质量守恒
		}

	
		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			//pMomentum += mParticles[i].mMass * velocity.norm();
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					double nodeMass = mGrid.getMass(grid_x, grid_y);
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx);
					mGrid.addVelocity(grid_x, grid_y, (particle.mMass * velocity + particle.mMass * particle.mCp * dx * mGrid.mCellSize) * weight / nodeMass);

				}
			}
		}


		for (int i = 0; i < mParticles.size(); i++)
		{
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			auto& Fe = particle.mFe;
			auto& Fp = particle.mFp;
			double Jp = particle.mDetFp;
			double Je = particle.mDetFe;
			double J = Jp * Je;
			auto FeT = Fe.transpose();
			double mu = 0.0;
			double lambda = 0.0;
			if (particle.mMaterials.mType == SNOW) {
				mu = particle.mMaterials.mu0 * std::exp(hardening * (1 - Jp));
				lambda = particle.mMaterials.lambda0 * std::exp(hardening * (1 - Jp));
			}
			else if(particle.mMaterials.mType == ELASITC){
				mu = particle.mMaterials.mu0;
				lambda = particle.mMaterials.lambda0;
			}

			JacobiSVD<Eigen::Matrix2d> svd(Fe,
				Eigen::ComputeFullU |
				Eigen::ComputeFullV);
			const Matrix2d& U = svd.matrixU();
			const Matrix2d& V = svd.matrixV();
			const Matrix2d& S = svd.singularValues().asDiagonal();
			const Matrix2d Re = U * V.transpose(); // 极分解 
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					// f = -Vp / J * (2*mu*(F-R)*FT + lambda * (J - 1) * J  * I) * dw;
					vec force = //  -particle.mVol * lambda * (Je - 1.0) * Je * dmat3(1.0) * dw;
						-particle.mVol * (2.0 * mu * (Fe - Re) * FeT + lambda * (Je - 1.0) * Je * Matrix2d::Identity()) * dw;
					mGrid.addForce(grid_x, grid_y, force);
				}
			}

		}


		for (int i = 0; i < mGrid.mNode.size(); i++) {
			mGrid.mNode[i].mForce += mGrid.mNode[i].mMass * mGravity;
			mGrid.mNode[i].mOldVelocity = mGrid.mNode[i].mVelocity;
			if (mGrid.mNode[i].mMass == 0.0) {
				mGrid.mNode[i].mVelocity(0) = 0.0;
				mGrid.mNode[i].mVelocity(1) = 0.0;
			}
			else {
				mGrid.mNode[i].mVelocity += mGrid.mNode[i].mForce / mGrid.mNode[i].mMass * dt;
			}
		}

		for (int i = 0; i < grid_num + 1; i++) {
			mGrid.mNode[i].mVelocity(1) = 0.0;
		}

		for (int i = 0; i < mParticles.size(); i++) {
			vec& pos = mParticles[i].mPosition;
			auto& particle = mParticles[i];
			vec& velocity = particle.mVelocity;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			Matrix2d gradV = Matrix2d::Zero();
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					vec dw = computeDerivativeWeightQuad(dx, mGrid.mCellSize);
					vec vel = mGrid.getVelocity(grid_x, grid_y);
					gradV += vel * dw.transpose();
				}
			}
			if (particle.mMaterials.mType == SNOW) {
				Matrix2d Fe_head = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
				Matrix2d F_next = Fe_head * particle.mFp;

				// clamp Fe_head
				Eigen::JacobiSVD<Eigen::Matrix2d> svd(
					Fe_head, Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::Vector2d sigma = svd.singularValues();
				double a = 1 - mCompression;
				double b = 1 + mStrech;
				sigma(0) = std::clamp(sigma(0), a, b);
				sigma(1) = std::clamp(sigma(1), a, b);
				particle.mFe = svd.matrixU() * sigma.asDiagonal() * svd.matrixV().transpose();
				particle.mFp = svd.matrixV() * sigma.cwiseInverse().asDiagonal() * svd.matrixU().transpose() * F_next;
				particle.mDetFe = particle.mFe.determinant();
				particle.mDetFp = particle.mFp.determinant();
			}
			else if (particle.mMaterials.mType == ELASITC) {
				particle.mFe = (Matrix2d::Identity() + dt * gradV) * particle.mFe;
				particle.mDetFe = particle.mFe.determinant();
			
			}

		}

		for (int i = 0; i < mParticles.size(); i++) {
			auto& particle = mParticles[i];
			vec V_apic = vec(0.0, 0.0);
			particle.mCp = Matrix2d::Zero();
			vec& pos = particle.mPosition;
			int c_x = static_cast<int>(floor((pos(0) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			int c_y = static_cast<int>(floor((pos(1) + 1.0) * 0.5 / mGrid.mCellSize - 0.5));
			vec posGrid = 0.5 * (pos + vec(1.0, 1.0));
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					int grid_x = c_x + r; // grid_x < 0;
					int grid_y = c_y + c;
					vec dx = (posGrid - mGrid.getGridPosition(grid_x, grid_y)) / mGrid.mCellSize;
					double weight = computeWeightQuad(dx); // 
					Node& node = mGrid.mNode[grid_y * (grid_num + 1) + grid_x];
					V_apic += weight * node.mVelocity;
					particle.mCp += 4 * weight / mGrid.mCellSize * node.mVelocity * dx.transpose();
				}
			}
			particle.mVelocity = V_apic;
		}

		for (auto& particle : mParticles) {
			vec p = particle.mPosition + dt * particle.mVelocity;
			if (p(1) < -0.8 || p(1) > 0.8) {
				particle.mVelocity(1) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
			if (p(0) < -0.8 || p(0) > 0.8) {
				particle.mVelocity(0) = 0.0; // particle.mVelocity - vy * n; // 把这个分量速度减掉
			}
		}

		for (auto& particle : mParticles) {
			particle.mPosition += dt * particle.mVelocity;
		}
		mGrid.clear();
	}

	void render_interaction(int frameIdx) {
		std::vector<Vector2d> pixels;
		std::vector<uint32_t> materials;
		std::vector<Vector3d> color;
		for (int i = 0; i < mParticles.size(); i++) {
			Vector2d ndcPos = mParticles[i].mPosition; //position 一定是 -1,1 之间的数
			pixels.push_back(ndcPos);
			materials.push_back(mParticles[i].mMaterials.mType);
		}

		mFrame.clear();
		mFrame.resize();
		// render ndc position
		for (int i = 0; i < pixels.size(); i++) {
			int pixel_x = static_cast<int>((pixels[i](0) + 1.0) * 0.5 * (mFrame.mWidth));
			int pixel_y = static_cast<int>(((pixels[i](1) + 1.0) * 0.5) * (mFrame.mHeight));
			if (materials[i] == SNOW) {
				for (int r = -1; r < 2; r++) {
					for (int c = -1; c < 2; c++) {
						mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
						mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 1);
						mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 1);
					}
				}
			}
			else if(materials[i] == ELASITC) {
				for (int r = -1; r < 2; r++) {
					for (int c = -1; c < 2; c++) {
						mFrame.drawColor(pixel_x + r, pixel_y + c, 0, 1);
						mFrame.drawColor(pixel_x + r, pixel_y + c, 1, 0);
						mFrame.drawColor(pixel_x + r, pixel_y + c, 2, 0);
					}
				}
			}

		}
		char buff[100];
		snprintf(buff, sizeof(buff), "interaction%03d.jpg", frameIdx);
		std::string filename = buff;
		mFrame.outputImg(filename);
	
	}

private:
	static double N(double dx) {
		double dx_abs = abs(dx);
		if (dx_abs < 0.5) {
			return 0.75 - dx_abs * dx_abs;
		}
		if (dx_abs < 1.5) {
			return 0.5 * (1.5 - dx_abs) * (1.5 - dx_abs);
		}
		return 0.0;
	}

	static double dN(double dx) {
		if (abs(dx) > 1.5) {
			return 0.0;
		}
		if (abs(dx) <= 0.5) {
			return -2.0 * dx;
		}
		else if(dx >= 0.5 && dx <= 1.5){
			return dx - 1.5;
		}
		else {
			return dx + 1.5;
		}
	}

	static double computeWeightQuad(vec dx) {
		return N(dx(0)) * N(dx(1));
	}

	static vec computeDerivativeWeightQuad(vec dx, double h) {
		return vec(dN(dx(0)) * N(dx(1)), N(dx(0)) * dN(dx(1))) / h;
	}

};