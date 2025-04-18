#include "common.h"

class Node {
public:
	Vector2d mVelocity;
	Vector2d mOldVelocity;
	double mMass;
	Vector2d mForce;

	Node() {
		mVelocity(0) = 0.0;
		mVelocity(1) = 0.0;
		mOldVelocity(0) = 0.0;
		mOldVelocity(1) = 0.0;
	}
	~Node() {
	}

};


class Grid2D {
public:
	uint32_t grid_num;
	double mCellSize;
	Grid2D() {}
	Grid2D(uint32_t ngrid) {
		grid_num = ngrid;
		mCellSize = 1.0 / grid_num;
		mNode.resize((grid_num + 1) * (grid_num + 1));
		for (int i = 0; i < mNode.size(); i++) {
			mNode[i].mForce(0) = 0.0;
			mNode[i].mForce(1) = 0.0;
			mNode[i].mMass = 0.0;
			mNode[i].mVelocity(0) = 0.0;
			mNode[i].mVelocity(1) = 0.0;
			mNode[i].mOldVelocity(0) = 0.0;
			mNode[i].mOldVelocity(1) = 0.0;
		}
	}


	Vector2d getGridPosition(int x, int y) {
		return Vector2d(x * mCellSize, y * mCellSize);
	}

	Vector2d& getVelocity(int x, int y) {
		return mNode[y * (grid_num + 1) + x].mVelocity;
	}

	void addMass(int x, int y, double mass) {
		mNode[y * (grid_num + 1) + x].mMass += mass;
	}

	double getMass(int x, int y) {
		return mNode[y * (grid_num + 1) + x].mMass;
	}

	void addVelocity(int x, int y, Vector2d vel) {
		mNode[y * (grid_num + 1) + x].mVelocity(0) += vel(0);
		mNode[y * (grid_num + 1) + x].mVelocity(1) += vel(1);

	}

	void clear() {
		for (int i = 0; i < mNode.size(); i++) {
			mNode[i].mForce(0) = 0.0;
			mNode[i].mForce(1) = 0.0;
			mNode[i].mMass = 0.0;
			mNode[i].mVelocity(0) = 0.0;
			mNode[i].mVelocity(1) = 0.0;
			mNode[i].mOldVelocity(0) = 0.0;
			mNode[i].mOldVelocity(1) = 0.0;
		}
	}

	void addForce(int x, int y,Vector2d& force) {
		mNode[y * (grid_num + 1) + x].mForce(0) += force(0);
		mNode[y * (grid_num + 1) + x].mForce(1) += force(1);
	}

	vector<Node> mNode;
private:
	

};