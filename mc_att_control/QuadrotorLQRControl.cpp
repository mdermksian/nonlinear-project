#include "QuadrotorLQRControl.hpp"

#include <conversion/rotation.h>
#include <drivers/drv_hrt.h>
#include <lib/ecl/geo/geo.h>
#include <circuit_breaker/circuit_breaker.h>
#include <mathlib/math/Limits.hpp>
#include <mathlib/math/Functions.hpp>

//#include <eigen3/Eigen/Dense> // ++ Math operations with Matrices ++

/* -------------------- added part ---------------------*/
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <memory>
#include <sstream>
#include "generate_reference.hpp"

using namespace matrix;
using namespace std;

string localpath("/cygdrive/c/PX4/home"); //<--- This is the path to the folder containing "Firmware"
// string localpath("/home/ubuntu/src"); //<--- This is the path to the folder containing "Firmware"

QuadrotorLQRControl::QuadrotorLQRControl()
{
    for (int i=0;i<nState;i++)
    {
        _current_state(i,0) = 0.0f;
        _current_state_ekf(i,0) = 0.0f;
        _eq_point(i,0) = 0.0f;
    }

    _ready_to_track = false;

    _eq_point(0,0) =  0.0f;
    _eq_point(1,0) =  0.0f;
    _eq_point(2,0) =  1.0f;

    _ref = _eq_point;

    ref_type = 0;

    u_control(0,0) = 0.0f;
    u_control(1,0) = 0.0f;
    u_control(2,0) = 0.0f;
    u_control(3,0) = 0.0f;

    gs_switch.loadRegions((localpath + "/Firmware/src/modules/mc_att_control/lqr_files/regions_switch.txt").c_str());
    gs_switch.loadControllers((localpath + "/Firmware/src/modules/mc_att_control/lqr_files/k_switch.txt").c_str());
    gs_switch.initializeRegion(_current_state(8, 0));

    gs_lin.loadRegions((localpath + "/Firmware/src/modules/mc_att_control/lqr_files/regions_lin.txt").c_str());
    gs_lin.loadControllers((localpath + "/Firmware/src/modules/mc_att_control/lqr_files/k_lin.txt").c_str());

    gs_contin.loadMatrices((localpath + "/Firmware/src/modules/mc_att_control/lqr_files/contin_param.txt").c_str());

    // _K = readMatrixK("/cygdrive/c/PX4/home/Firmware/src/modules/mc_att_control/lqr_files/new_controller.txt");
    // _PMATRIX = readMatrixP("/cygdrive/c/PX4/home/Firmware/src/modules/mc_att_control/lqr_files/new_pe.txt");

    // ff_thrust = 7.848f; //  [N]
    ff_thrust = 0.8f * 9.81f;   //  [N]

    _auto_eq_point_flag = true;

    ofstream outfile1;
    outfile1.open((localpath + "/Firmware/src/modules/mc_att_control/output_files/control_input.txt").c_str(), std::ios::out);
    outfile1.close();

    ofstream outfile3;
    outfile3.open((localpath + "/Firmware/src/modules/mc_att_control/output_files/state.txt").c_str(), std::ios::out);
    outfile3.close();

    ofstream outfile5;
    outfile5.open((localpath + "/Firmware/src/modules/mc_att_control/output_files/lyapunov.txt").c_str(), std::ios::out);
    outfile5.close();

    ofstream outfile4;
    outfile4.open((localpath + "/Firmware/src/modules/mc_att_control/output_files/ekf.txt").c_str(), std::ios::out);
    outfile4.close();

    ofstream outfile2;
    outfile2.open((localpath + "/Firmware/src/modules/mc_att_control/output_files/ref.txt").c_str(), std::ios::out);
    outfile2.close();

    _past_time = hrt_absolute_time() * 1e-6;
}



Matrix <float, nCont, nState>  QuadrotorLQRControl::readMatrixK(const char *filename)
{
    cout << "TRYING TO READ K\n";
    static Matrix <float, nCont, nState> result;
    static int rows = nCont;
    static int cols = nState;
    ifstream infile;
    infile.open(filename);
    if (infile.is_open()){
        for (int i=0; i<rows;i++){
    		string line;
    		getline(infile, line);
    		stringstream stream(line);
    		for (int j=0; j<cols; j++){
    			stream >> result(i,j);
    		}

    	}
    	infile.close();
    }else cout << "Unable to open file\n";
    return result;
}

bool QuadrotorLQRControl::isReadyForTracking(void)
{
    return _ready_to_track;
}

Matrix<float,nCont,1> QuadrotorLQRControl::LQRcontrol()
{
    // Select state origin (true or EKF)
    // Matrix<float,nState,1> state = _current_state;  
    Matrix<float,nState,1> state = _current_state_ekf;

    if(!_ready_to_track){
        // Check if within 1cm of equilibrium point euclidean distance (not necessarily stable!)
        float dist = (pow(state(0,0)-_eq_point(0,0), 2) + pow(state(1,0)-_eq_point(1,0), 2) + pow(state(2,0)-_eq_point(2,0),2));
        // if( dist < 1e-3 && fabs(state(14, 0)) < 0.05f) {
        if( dist < 1e-3 && fabs(state(6, 0)) < 0.005f) {
            _ready_to_track = true;
        }
    }

    // cout << _ready_to_track << endl;

    //static Matrix<float,4,1> u_control;
    static Matrix<float,nCont,1> u_control_norm;
    static Matrix<float,nState,1> delta_x;
    static Matrix<float, 1,nState> v_b;
    static Matrix<float,1,nState> delta_x_T;
    static Matrix<float,1,1> _lyap_fun;     
    const hrt_abstime now = hrt_absolute_time();

    float _current_time = now *1e-6;
    // float dt = _current_time-_past_time;
     
    _past_time = _current_time;

    if(_ready_to_track){
        // delta_x = _current_state - _ref;
        // cout << _ref(0, 0) << ", " << _ref(1, 0) << ", " << _ref(2, 0) << endl;
        delta_x = state - _ref;
    } else {
        // delta_x   = _current_state - _eq_point;  
        delta_x   = state - _eq_point;   
    }

    // cout << state(12,0) << ", " << state(13,0) << ", " << state(14,0) << ", " << state(15,0) << endl;

    // Matrix<float,nCont,nState> K = gs_switch.getK(state(8,0));
    // Matrix<float,nCont,nState> K = gs_lin.getK(state(8,0));
    Matrix<float,nCont,nState> K = gs_contin.getK(state(8,0));
    u_control = -K*(delta_x);
    
    // cout << gs_lin.getRegionInd(state(8,0)) << ", " << state(8,0) << endl;

    // delta_x_T = delta_x.transpose();
    
    // v_b = delta_x_T*_PMATRIX;
    // _lyap_fun = v_b*delta_x;
    //cout<< dt << "\t" << _P(0,0) << "\n";

    float phi = state(6,0);
    float theta = state(7,0);
    float ff_control = ff_thrust * cos(phi) * cos(theta); // scale gravity feedforward based on overall tilt
   // !! IMPORTANT scale the control inputs.......


    // u_control_norm(0,0) = fmin(fmax((u_control(0,0)+ff_thrust)/16.0f, 0.0f), 1.0f);
    u_control_norm(0,0) = fmin(fmax((u_control(0,0)+ff_control)/16.0f, 0.0f), 1.0f);
    u_control_norm(1,0) = fmin(fmax((u_control(1,0))/(4.0f), -1.0f), 1.0f);  
    u_control_norm(2,0) = fmin(fmax((u_control(2,0))/(4.0f), -1.0f), 1.0f);
    // u_control_norm(3,0) = fmin(fmax((u_control(3,0))/(1.0f), -1.0f), 1.0f);
    u_control_norm(3,0) = fmin(fmax((u_control(3,0))/(0.05f), -1.0f), 1.0f);

    // cout << u_control_norm(0,0) << ", " << u_control_norm(1,0) << ", " << u_control_norm(2,0) << ", " << u_control_norm(3,0) << endl; 

   // not normalized control inputs
     u_control(0,0) = u_control_norm(0,0)*16.0f;
     u_control(1,0) = u_control_norm(1,0)*4.0f;
     u_control(2,0) = u_control_norm(2,0)*4.0f;
     u_control(3,0) = u_control_norm(3,0)*0.05f;
     
    //"\t" <<  u_control(0,0)+ff_thrust << "\n";
         /* Save data*/
    if(_ready_to_track) {
        writeStateOnFile((localpath + "/Firmware/src/modules/mc_att_control/output_files/state.txt").c_str(), _current_state, now);
        writeInputOnFile((localpath + "/Firmware/src/modules/mc_att_control/output_files/control_input.txt").c_str(), u_control_norm, now); 
        writeLyapunovOnFile((localpath + "/Firmware/src/modules/mc_att_control/output_files/lyapunov.txt").c_str(), _lyap_fun(0,0), now); 
        writeStateOnFile((localpath + "/Firmware/src/modules/mc_att_control/output_files/ekf.txt").c_str(), _current_state_ekf, now);
        writeReferenceOnFile((localpath + "/Firmware/src/modules/mc_att_control/output_files/ref.txt").c_str(), _ref, now);
    }
    
    return u_control_norm;    

}

Matrix<float,nCont,1> QuadrotorLQRControl::normalizationControlInputs(Matrix<float,nCont,1> _u)
{
   Matrix<float,nCont,1> _u_norm;
   _u_norm(0,0) = _u(0,0)*16.0f;
   _u_norm(1,0) = _u(1,0)*(0.1080f*4.0f);
   _u_norm(2,0) = _u(2,0)*(0.1080f*4.0f);
   _u_norm(3,0) = _u(3,0)*(0.1f*1.0f); 

    return _u_norm;
}

Matrix<float,nCont,1> QuadrotorLQRControl::getLQRcontrols()
{
    return u_control;
}

void QuadrotorLQRControl::setCurrentState(Matrix<float,nState,1> _state_estimate)
{
    _current_state = _state_estimate;
}

void QuadrotorLQRControl::computeIntegral(Matrix<float, nState, 1> &state)
{
    const hrt_abstime now = hrt_absolute_time();
    float _current_time = now *1e-6;
    float dt = _current_time - _past_time;

    if(dt > 1.0f) return;

    int integ_states[] = {0, 1, 2, 8};
    
    Matrix<float,4,1> sigma;
    for(size_t i = 0; i < nRef; ++i) {
        int ind = integ_states[i];
        float error;
        if(_ready_to_track) {
            error = state(ind,0) - _ref(ind,0);
        } else {
            error = state(ind,0) - _eq_point(ind,0);
        }
        
        float sigma0 = state(12+i,0);
        float dsigma0 = error;
        float sigma_mid = sigma0 + dt * dsigma0;
        state(12+i,0) = sigma_mid;
        // Matrix<float,4,1> refs_i;
        // if(_ready_to_track){
        //     refs_i = generate_reference(_current_time+dt, ref_type, _ref(8,0));
        // } else {
        //     refs_i = _eq_point.slice<4,1>(0,0);
        //     refs_i(3, 0) = _eq_point(8, 0);
        // }
        // float ref_i = refs_i(i,0);
        // float dsigma_mid = sigma_mid - ref_i;
        // state(12+i,0) = sigma0 + dt * (dsigma0 + dsigma_mid) / 2;
    }
}

void QuadrotorLQRControl::setCurrentState(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos)
{
    // _current_state(0,0) = _v_local_pos.x;
    // _current_state(1,0) = _v_local_pos.y;
    // _current_state(2,0) = -_v_local_pos.z;
    // _current_state(3,0) = _v_local_pos.vx;
    // _current_state(4,0) = _v_local_pos.vy;
    // _current_state(5,0) = _v_local_pos.vz;

    // _current_state(6,0)  = Eulerf(Quatf(_v_att.q)).phi();
    // _current_state(7,0)  = Eulerf(Quatf(_v_att.q)).theta();
    // _current_state(8,0) = Eulerf(Quatf(_v_att.q)).psi();
    // _current_state(9,0)  = _v_att.rollspeed;
    // _current_state(10,0)  = _v_att.pitchspeed;
    // _current_state(11,0) = _v_att.yawspeed;

    _current_state(0,0) = _v_local_pos.x;
    _current_state(1,0) = _v_local_pos.y;
    _current_state(2,0) = -1 * _v_local_pos.z;
    _current_state(6,0) = Eulerf(Quatf(_v_att.q)).phi();
    _current_state(7,0) = Eulerf(Quatf(_v_att.q)).theta();
    float psi = Eulerf(Quatf(_v_att.q)).psi();
    float last_psi = _current_state(8,0);
    psi = unwrap2pi(psi, last_psi);
    _current_state(8,0) = psi;
    _current_state(9,0) = _v_att.rollspeed;
    _current_state(10,0) = _v_att.pitchspeed;
    _current_state(11,0) = _v_att.yawspeed;

    Matrix<float,3,1> world_vel;
    world_vel(0,0) = _v_local_pos.vx;
    world_vel(1,0) = _v_local_pos.vy;
    world_vel(2,0) = _v_local_pos.vz;

    Matrix<float,3,1> euler_ang = _current_state.slice<3,1>(6,0);
    Matrix<float,3,3> rot_matrix = world_to_body_rot(euler_ang);

    Matrix<float,3,1> body_vel = rot_matrix * world_vel;
    _current_state(3,0) = body_vel(0,0);
    _current_state(4,0) = body_vel(1,0);
    _current_state(5,0) = body_vel(2,0);

    computeIntegral(_current_state);
}

void QuadrotorLQRControl::setCurrentStateEkf(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos)
{
    // _current_state_ekf(0,0) = _v_local_pos.x;
    // _current_state_ekf(1,0) = _v_local_pos.y;
    // _current_state_ekf(2,0) = -_v_local_pos.z;
    // _current_state_ekf(3,0) = _v_local_pos.vx;
    // _current_state_ekf(4,0) = _v_local_pos.vy;
    // _current_state_ekf(5,0) = _v_local_pos.vz;

    // _current_state_ekf(6,0)  = Eulerf(Quatf(_v_att.q)).phi();
    // _current_state_ekf(7,0)  = Eulerf(Quatf(_v_att.q)).theta();
    // _current_state_ekf(8,0) = Eulerf(Quatf(_v_att.q)).psi();
    // _current_state_ekf(9,0)  = _v_att.rollspeed;
    // _current_state_ekf(10,0)  = _v_att.pitchspeed;
    // _current_state_ekf(11,0) = _v_att.yawspeed;

    _current_state_ekf(0,0) = _v_local_pos.x;
    _current_state_ekf(1,0) = _v_local_pos.y;
    _current_state_ekf(2,0) = -1 * _v_local_pos.z;
    _current_state_ekf(6,0) = Eulerf(Quatf(_v_att.q)).phi();
    _current_state_ekf(7,0) = Eulerf(Quatf(_v_att.q)).theta();
    float psi = Eulerf(Quatf(_v_att.q)).psi();
    float last_psi = _current_state_ekf(8,0);
    psi = unwrap2pi(psi, last_psi);
    _current_state_ekf(8,0) = psi;
    _current_state_ekf(9,0) = _v_att.rollspeed;
    _current_state_ekf(10,0) = _v_att.pitchspeed;
    _current_state_ekf(11,0) = _v_att.yawspeed;

    Matrix<float,3,1> world_vel;
    world_vel(0,0) = _v_local_pos.vx;
    world_vel(1,0) = _v_local_pos.vy;
    world_vel(2,0) = _v_local_pos.vz;

    Matrix<float,3,1> euler_ang = _current_state_ekf.slice<3,1>(6,0);
    Matrix<float,3,3> rot_matrix = world_to_body_rot(euler_ang);

    Matrix<float,3,1> body_vel = rot_matrix * world_vel;
    _current_state_ekf(3,0) = body_vel(0,0);
    _current_state_ekf(4,0) = body_vel(1,0);
    _current_state_ekf(5,0) = body_vel(2,0);
    
    computeIntegral(_current_state_ekf);
}


void QuadrotorLQRControl::setAutoEqPoint(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos)
{
    _eq_point(0,0) = _v_local_pos.x;
    _eq_point(1,0) = _v_local_pos.y;
    _eq_point(2,0) = -_v_local_pos.z;
    _eq_point(3,0) = _v_local_pos.vx;
    _eq_point(4,0) = _v_local_pos.vy;
    _eq_point(5,0) = _v_local_pos.vz;

    _eq_point(6,0)  = Eulerf(Quatf(_v_att.q)).phi();
    _eq_point(7,0)  = Eulerf(Quatf(_v_att.q)).theta();
    _eq_point(8,0) = Eulerf(Quatf(_v_att.q)).psi();
    _eq_point(9,0)  = _v_att.rollspeed;
    _eq_point(10,0)  = _v_att.pitchspeed;
    _eq_point(11,0) = _v_att.yawspeed;	
}

void QuadrotorLQRControl::setEquilibriumPoint(Matrix<float,nState,1> eqPoint)
{
    _eq_point(0,0) = eqPoint(0,0);
    _eq_point(1,0) = eqPoint(1,0);
    _eq_point(2,0) = eqPoint(2,0);
    _eq_point(3,0) = eqPoint(3,0);
    _eq_point(4,0) = eqPoint(4,0);
    _eq_point(5,0) = eqPoint(5,0);

    _eq_point(6,0)  = eqPoint(6,0);
    _eq_point(7,0)  = eqPoint(7,0);
    _eq_point(8,0)  = eqPoint(8,0);
    _eq_point(9,0)  = eqPoint(9,0);
    _eq_point(10,0) = eqPoint(10,0);
    _eq_point(11,0) = eqPoint(11,0);
}

void QuadrotorLQRControl::setReferencePoint(Matrix<float,4,1> ref)
{
    _ref(0, 0) = ref(0, 0);
    _ref(1, 0) = ref(1, 0);
    _ref(2, 0) = ref(2, 0);
    _ref(8, 0) = ref(3, 0);
}

void QuadrotorLQRControl::setReferenceType(int type) {
    ref_type = type;
}

void QuadrotorLQRControl::setAutoEqPointFlag(bool flag)
{
    _auto_eq_point_flag = flag;
}

bool QuadrotorLQRControl::getAutoEqPointFlag()
{
    return _auto_eq_point_flag;
}

/* Save data on files */

void QuadrotorLQRControl::writeStateOnFile(const char *filename, Matrix <float, nState, 1> vect, hrt_abstime t) {

	ofstream outfile;
	outfile.open(filename, std::ios::out | std::ios::app);
    outfile << t << "\t";   // time
    
	for(int i=0;i<nState;i++){
		if(i==nState-1){
			outfile << vect(i,0) << "\n";
		}else{
            outfile << vect(i,0) << "\t";
		}
	}
	outfile.close();
	return;
}


void QuadrotorLQRControl::writeInputOnFile(const char *filename, Matrix <float, nCont, 1> vect, hrt_abstime t) {

	ofstream outfile;
	outfile.open(filename, std::ios::out | std::ios::app);
    outfile << t << "\t";   // time
        
	for(int i=0;i<nCont;i++){
		if(i==3){
			outfile << vect(i,0) << "\n";
		}else{
            outfile << vect(i,0) << "\t";
		}
	}
	outfile.close();
	return;
}

void QuadrotorLQRControl::writeLyapunovOnFile(const char *filename, float value, hrt_abstime t) {

	ofstream outfile;
	outfile.open(filename, std::ios::out | std::ios::app);
    outfile << t << "\t" << value << "\n";   
	outfile.close();
	return;
}

void QuadrotorLQRControl::writeReferenceOnFile(const char *filename, Matrix <float, nState, 1> vect, hrt_abstime t) {
    const int inds[] = {0, 1, 2, 8};

    ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::app);
    outfile << t << "\t";   // time

        
    for(int i=0;i<nRef;i++){
        if(i==nRef-1){
            outfile << vect(inds[i],0) << "\n";
        }else{
            outfile << vect(inds[i],0) << "\t";
        }
    }
    outfile.close();
    return;
}

Matrix <float, nState, nState> QuadrotorLQRControl::readMatrixP(const char *filename)
{
    cout << "TRYING TO READ P\n";
    static Matrix <float, nState, nState> result;
    static int rows = nState;
    static int cols = nState;
    ifstream infile;
    infile.open(filename);
    if (infile.is_open()){
        for (int i = 0; i<rows; i++){
            string line;
            getline(infile, line);
            stringstream stream(line);
            for (int j=0; j<cols; j++){
                stream >> result(i,j);
            }
        }
        infile.close();
    }else cout << "Unable to open file\n";
    return result;
}

// Matrix <float, 12, 12>  QuadrotorLQRControl::readMatrixP(const char *filename)
// {
//     static Matrix <float, 12, 12> result;
//     static int rows = 12;
//     static int cols = 12;
//     ifstream infile;
//     infile.open(filename);
//     if (infile.is_open()){
//          for (int i=0; i<rows;i++){
//     		string line;
//     		getline(infile, line);
//     		stringstream stream(line);
//     		for (int j=0; j<cols; j++){
//     			stream >> result(i,j);
//     		}

//     	}
//     	infile.close();
//     }else cout << "Unable to open file";
//     return result;

//  }

Matrix<float,3,3> QuadrotorLQRControl::world_to_body_rot(Matrix<float,3,1> euler_ang) {
  float phi = euler_ang(0,0);
  float theta = euler_ang(1,0);
  float psi = euler_ang(2,0);

  Matrix<float,3,3> rot_mat;

  rot_mat(0,0) = cos(theta) * cos(psi);
  rot_mat(0,1) = cos(theta) * sin(psi);
  rot_mat(0,2) = -1 * sin(theta);
  rot_mat(1,0) = sin(phi) * sin(theta) * cos(psi)  -  cos(phi) * sin(psi);
  rot_mat(1,1) = sin(phi) * sin(theta) * sin(psi)  +  cos(phi) * cos(psi);
  rot_mat(1,2) = sin(phi) * cos(theta);
  rot_mat(2,0) = cos(phi) * sin(theta) * cos(psi)  +  sin(phi) * sin(psi);
  rot_mat(2,1) = cos(phi) * sin(theta) * sin(psi)  -  sin(phi) * cos(psi);
  rot_mat(2,2) = cos(phi) * cos(theta);

  return rot_mat;
}
