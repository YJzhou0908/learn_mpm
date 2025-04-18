#include "common.h"

class Materials {
public:
	Materials() {}
	Materials(uint32_t Type){
		if (Type == ELASITC) {
			mE = 1.4e5;
			mNu = 0.2;
			hardening = 2.0;
			flipAlpha = 0.95;
			mCompression = 5.5e-2;
			mStrech = 7.5e-3;
			mu0 = mE / (2 * (1 + mNu));
			lambda0 = (mE * mNu) / ((1 + mNu) * (1 - 2 * mNu));
			mType = ELASITC;
		}
		else if (Type == SNOW) {
			mE = 1.4e5;
			mNu = 0.2;
			hardening = 2.0;
			flipAlpha = 0.95;
			mCompression = 5.5e-2;
			mStrech = 7.5e-3;
			mu0 = mE / (2 * (1 + mNu));
			lambda0 = (mE * mNu) / ((1 + mNu) * (1 - 2 * mNu));
			mType = SNOW;
		
		}
	}

	double mE ;
	double mNu ; // Possion ratio
	double hardening; // 2.0 ºÍ 7.5e-2
	double flipAlpha ;
	double ground;
	double mCompression;
	double mStrech;
	double mu0;
	double lambda0;

	uint32_t mType;
};