#include <cmath>
#include <matrix/matrix/math.hpp>
#include "generate_reference.hpp"

using namespace matrix;

Matrix<float,4,1>  generate_reference(float t, int type, float last_psir){
	Matrix<float,4,1> pt;
	float PI = 3.14159265;
	if(type == 0){ // Helix
		// float rad = 0.25f;
		float rad = 0.0f;
		float w = 0.1f;
		// float z_rate = 0.5f;
		float z_rate = 0.0f;

		float psir = 2*PI * w * t;
		float pNr = rad * cos(psir);
		float pEr = rad * sin(psir);
		float hr = z_rate * t;

		pt(0,0) = pNr;
		pt(1,0) = pEr;
		pt(2,0) = hr+1.0f;
		pt(3,0) = psir;
	} else { // Lissajous
		float A = 2.0f;
		float a = 0.1f;
		float d = 0.0f;
		float B = 2.0f;
		float b = 0.15f;
		float z_step = 1.0f;
		float z_per = 10.0f;

		float dpNr = 2*PI * a * A * cos(2*PI*a*t + d);
    	float dpEr = 2*PI * b * B * cos(2*PI*b*t);
		float psir = atan2(dpEr, dpNr);
		psir = unwrap2pi(psir, last_psir);
		float pNr = A * sin(2*PI*a*t + d);
    	float pEr = B * sin(2*PI*b*t);
		float hr = z_step * floor(t / z_per);

		pt(0,0) = pNr;
		pt(1,0) = pEr;
		pt(2,0) = hr+1.0f;
		pt(3,0) = psir;
	}

	return pt;
};

float unwrap2pi(float ang, float lastang) {
	const int PI = 3.14159265;
	while(abs(ang - lastang) > PI) {
		float val = ang - lastang;
		int sign = (float(0) < val) - (val < float(0));
		ang = ang - 2*PI * sign;
	}
	return ang;
}