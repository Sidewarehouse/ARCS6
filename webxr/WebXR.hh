#ifndef WEB_XR
#define WEB_XR

#include "ConstParams.hh"
#include "WebXRMotorModel.hh"

#define ACTUATOR_STATUS_NORMAL 1
#define ACTUATOR_STATUS_INACTIVE 2

extern "C" void setActuatorStatus(int axis, int status, double current,
                                  double position);
extern "C" void setReadVariable(int number, double value);
extern "C" bool canGetWriteVariable(int number);
extern "C" double getWriteVariable(int number);
extern "C" void setWriteVariable(int number, double value);
extern "C" void setPlot(int figureNumber, int variableNumber, double value);
extern "C" void setAxisRadian(int axisNumber, double radian);
extern "C" double getAxisRadian(int axisNumber);
extern "C" void debug(double value);

extern long currentTimeInMsec;
extern double currentTime;
extern MotorModel MotorModels[ARCS::ConstParams::ACTUATOR_NUM];

#endif
