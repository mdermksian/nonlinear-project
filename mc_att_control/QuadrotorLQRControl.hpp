#include <lib/mixer/mixer.h>
#include <mathlib/math/filter/LowPassFilter2pVector3f.hpp>
#include <matrix/matrix/math.hpp>
#include <perf/perf_counter.h>
#include <px4_config.h>
#include <px4_defines.h>
#include <px4_module.h>
#include <px4_module_params.h>
#include <px4_posix.h>
#include <px4_tasks.h>
#include <uORB/topics/actuator_controls.h>
#include <uORB/topics/battery_status.h>
#include <uORB/topics/manual_control_setpoint.h>
#include <uORB/topics/multirotor_motor_limits.h>
#include <uORB/topics/parameter_update.h>
#include <uORB/topics/rate_ctrl_status.h>
#include <uORB/topics/sensor_bias.h>
#include <uORB/topics/sensor_correction.h>
#include <uORB/topics/sensor_gyro.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_attitude_setpoint.h>
#include <uORB/topics/vehicle_control_mode.h>
#include <uORB/topics/vehicle_rates_setpoint.h>
#include <uORB/topics/vehicle_status.h>
#include <uORB/topics/vehicle_land_detected.h>
#include <uORB/topics/landing_gear.h>
#include <drivers/drv_hrt.h>

#include <uORB/topics/experiment_mode.h>
#include <uORB/topics/vehicle_local_position.h>

#include "GainSchedule.hpp"

using namespace matrix;

static const int nState = 16;
static const int nCont = 4;
static const int nRef = 4;

class QuadrotorLQRControl
{
public:
       QuadrotorLQRControl();
       
       void setEquilibriumPoint(Matrix<float,nState,1> eqPoint);

       void setReferencePoint(Matrix<float,nRef, 1> ref);

       void setCurrentState(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos);
       
       void setCurrentStateEkf(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos);
       
       Matrix<float,nCont,1> LQRcontrol();
       
       void setAutoEqPointFlag(bool flag);
    
       void setAutoEqPoint(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos); 

       bool getAutoEqPointFlag();

       Matrix<float,nCont,1> getLQRcontrols();
       
       void setCurrentState(Matrix<float,nState,1> _state_estimate);

       Matrix<float,nCont,1> normalizationControlInputs(Matrix<float,nCont,1> _u);

       bool isReadyForTracking(void);

       void setReferenceType(int type);
       

private:

      bool _ready_to_track;
       
      GainScheduleLin gs_lin;

      GainScheduleSwitch gs_switch;

      GainScheduleContin gs_contin;
       
      Matrix<float,nState,1> _eq_point;

      Matrix<float,nState,1> _ref;

      int ref_type;

      Matrix<float,nCont,nState> _K;
      
      Matrix<float,nState,nState> _PMATRIX;
      
      Matrix<float,nState,1> _current_state;

      Matrix<float,nState,1> _current_state_ekf;


      bool _auto_eq_point_flag;
     
      float ff_thrust;

      // void computeIntegral(struct vehicle_attitude_s _v_att, struct vehicle_local_position_s  _v_local_pos);
      void computeIntegral();

      void writeStateOnFile(const char *filename, Matrix <float, nState, 1> vect, hrt_abstime t); 

      Matrix <float, nCont, nState> readMatrixK(const char *filename);

      void writeInputOnFile(const char *filename, Matrix <float, nCont, 1> vect, hrt_abstime t); 

      void writeLyapunovOnFile(const char *filename, float value, hrt_abstime t);

      Matrix <float, nState, nState> readMatrixP(const char *filename);

      float _past_time;

      Matrix<float,nCont,1> u_control;
      
      Matrix<float,3,3> world_to_body_rot(Matrix<float,3,1> euler_ang);
};
