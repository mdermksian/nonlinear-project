#include "GainSchedule.hpp"
#include <matrix/matrix/math.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

/********************************************/
/* LINEAR INTERPOLATION GAIN SCHEDULE CLASS */
/********************************************/
GainScheduleLin::GainScheduleLin() {
    min = 0.0;
    max = 0.0;
    numreg = 0;
}

void GainScheduleLin::loadRegions(const char *filename) {
    std::cout << "Loading regions\n";
    std::ifstream infile;
    infile.open(filename);
    if(infile.is_open()) {
        while(true){
            std::string line;
            std::getline(infile, line);
            float point = std::stof(line);
            regions.push_back(point);
            if(infile.eof()) break;
        }
    } else std::cout << "Unable to open regions file\n";

    min = regions.front();
    max = regions.back();
    numreg = regions.size()-1;
}

void GainScheduleLin::loadControllers(const char *filename) {
    std::cout << "Loading controllers\n";
    std::ifstream infile;
    infile.open(filename);
    if(infile.is_open()) {
        while(true){
            static matrix::Matrix <float,nRow,nCol> result;
            for(int i = 0; i < nRow; ++i){
                std::string line;
                std::getline(infile, line);
                std::stringstream stream(line);
                std::cout << line << std::endl;
                for(int j = 0; j < nCol; ++j){
                    stream >> result(i,j);
                }
                if(infile.eof()) break;
            }
            K.push_back(result);
            if(infile.eof()) break;
        }
    } else std::cout << "Unable to open controllers file\n";

    if(K.size() != regions.size()) std::cout << "WARNING!!!!! NUMBER OF CONTROLLERS DOES NOT MATCH NUMBER OF REGIONS\n";
}

int GainScheduleLin::getRegionInd(float psi) {
    float range = max - min;
    float step = range/(numreg);
    int ind = (psi-min)/step;
    if(ind >= numreg) {
        ind = numreg - 1;
    } else if(ind < 0) {
        ind = 0;
    }
    return ind;
}

matrix::Matrix<float,nRow,nCol> GainScheduleLin::getK(float psi) {
    int ind = getRegionInd(psi);
    float pos = (psi - regions[ind])/(regions[ind+1] - regions[ind]);
    matrix::Matrix<float,nRow,nCol> K1 = K[ind];
    matrix::Matrix<float,nRow,nCol> K2 = K[ind+1];

    return pos*K2+(1-pos)*K1;
}

/********************************************/
/*      SWITCHING GAIN SCHEDULE CLASS       */
/********************************************/
GainScheduleSwitch::GainScheduleSwitch() {
    min = 0.0;
    max = 0.0;
    numreg = 0;
    cur_reg = 0;
}

void GainScheduleSwitch::loadRegions(const char *filename) {
    std::cout << "Loading regions\n";
    std::ifstream infile;
    infile.open(filename);
    if(infile.is_open()) {
        while(true){
            std::string line;
            std::getline(infile, line);
            float temp[2];
            std::stringstream stream(line);
            for(int i = 0; i < 2; ++i){
                stream >> temp[i];
            }
            Region tempreg(temp);
            tempreg.print();
            regions.push_back(tempreg);
            if(infile.eof()) break;
        }
    } else std::cout << "Unable to open regions file\n";

    min = regions.front().start;
    max = regions.back().end;
    numreg = regions.size();
}

void GainScheduleSwitch::loadControllers(const char *filename) {
    std::cout << "Loading controllers\n";
    std::ifstream infile;
    infile.open(filename);
    if(infile.is_open()) {
        while(true){
            static matrix::Matrix <float,nRow,nCol> result;
            for(int i = 0; i < nRow; ++i){
                std::string line;
                std::getline(infile, line);
                std::stringstream stream(line);
                for(int j = 0; j < nCol; ++j){
                    stream >> result(i,j);
                }
                if(infile.eof()) break;
            }
            K.push_back(result);
            if(infile.eof()) break;
        }
    } else std::cout << "Unable to open controllers file\n";

    if(K.size() != regions.size()) std::cout << "WARNING!!!!! NUMBER OF CONTROLLERS DOES NOT MATCH NUMBER OF REGIONS\n";
}

void GainScheduleSwitch::initializeRegion(float psi) {
    int ind = -1;
    float lowest = 1e10; //Arbitrary large number to ensure that we find a lower one
    for(int i = 0; i < (int)regions.size(); ++i) {
        float center = (regions[i].end - regions[i].start)/2;
        float center_dist = abs(center-psi);
        if(center_dist < lowest) {
            lowest = center_dist;
            ind = i;
        }
    }

    if(ind < 0 || ind >= (int)regions.size()){
        std::cout << "Something went wrong initializing the region";
        cur_reg = 0;
    } else {
        cur_reg = ind;
    }
}

int GainScheduleSwitch::getRegionInd(float psi) {
    int ind = cur_reg;
    Region current = regions[cur_reg];
    if(psi > current.end) {
        ind = cur_reg+1;
    } else if(psi < current.start) {
        ind = cur_reg-1;
    }

    if(ind < 0) {
        ind = 0;
    } else if (ind >= (int)regions.size()){
        ind = (int)regions.size()-1;
    }

    cur_reg = ind;
    return ind;
}

matrix::Matrix<float,nRow,nCol> GainScheduleSwitch::getK(float psi) {
    int ind = getRegionInd(psi);
    return K[ind];
}

/********************************************/
/*      CONTINUOUS GAIN SCHEDULE CLASS      */
/********************************************/
GainScheduleContin::GainScheduleContin(){

}

matrix::Matrix<float,nRow,nCol> GainScheduleContin::getK(float x) {

}


/*******************************************/
/*              REGION CLASS               */
/*******************************************/
Region::Region() {
    start = 0.0;
    end = 0.0;
}

Region::Region(float s, float e) {
    start = s;
    end = e;
}

Region::Region(float pair[2]) {
    start = pair[0];
    end = pair[1];
}

void Region::print(void) {
    std::cout << "(" << start << ", " << end << ")\n";
}