#include <cmath>
#include <matrix/matrix/math.hpp>

using namespace matrix;

Matrix<float,4,1>  generate_reference(float t, int type, float last_psir);

float wrap2pi(float ang, float lastang);