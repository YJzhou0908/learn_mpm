#include "common.h"
#include "Materials.h"

class Particle2D {
public:
	Vector2d mPosition;
	Vector2d mVelocity;
	double mMass;
	double mVol;
	double mDensity;

	Materials mMaterials;

	Vector3d mColor;

	// deforamtion gradiant, elastic and plastic
	Matrix2d mFe, mFp;
	Matrix2d mCp; // APIC 
	double mDetFe, mDetFp;

	Particle2D() {}
	Particle2D(Vector2d& pos, Vector2d& vel, double mass, double vol, double density, uint32_t material_type) {
		mPosition = pos;
		mVelocity = vel;
		mMass = mass;
		mVol = vol;
		mDensity = density;
		mFe = Matrix2d::Identity();
		mFp = Matrix2d::Identity();
		mCp = Matrix2d::Zero();
		mDetFe = 1.0;
		mDetFp = 1.0;
		mMaterials = Materials(material_type);

	}
	Particle2D(Vector2d& pos, Vector2d& vel, double mass, double vol, double density) {
		mPosition = pos;
		mVelocity = vel;
		mMass = mass;
		mVol = vol;
		mDensity = density;
		mFe = Matrix2d::Identity();
		mFp = Matrix2d::Identity();
		mCp = Matrix2d::Zero();
		mDetFe = 1.0;
		mDetFp = 1.0;
	}
	~Particle2D() {
	}

private:


};