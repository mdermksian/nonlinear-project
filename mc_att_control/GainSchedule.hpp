#include <matrix/matrix/math.hpp>
#include <vector>

// Irresponsibly defined global constants
static const int nRow = 4;
static const int nCol = 16;

/* REGION UTILITY CLASS */
class Region {
public:
    float start;
    float end;
    Region();
    Region(float s, float e);
    Region(float pair[2]);
    void print(void);
};

/* LINEAR INTERPOLATED GAIN SCHEDULER */
class GainScheduleLin {
public:
    GainScheduleLin();
    void loadRegions(const char *filename);
    void loadControllers(const char *filename);
    int getRegionInd(float psi);
    matrix::Matrix<float,nRow,nCol> getK(float psi);

private:
    float min;
    float max;
    int numreg;
    std::vector<float> regions;
    std::vector<matrix::Matrix<float,nRow,nCol>> K;

};

/* SWITCHED MODE GAIN SCHEDULER */
class GainScheduleSwitch {
public:
    GainScheduleSwitch();
    void loadRegions(const char *filename);
    void loadControllers(const char *filename);
    void initializeRegion(float psi);
    int getRegionInd(float psi);
    matrix::Matrix<float,nRow,nCol> getK(float psi);

private:
    float min;
    float max;
    int numreg;
    int cur_reg;
    std::vector<Region> regions;
    std::vector<matrix::Matrix<float,nRow,nCol>> K;
};

/* CONTINUOUS MODE GAIN SCHEDULER */
class GainScheduleContin {
public:
    GainScheduleContin();
    matrix::Matrix<float,nRow,nCol> getK(float x);
};