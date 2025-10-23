#include "WebXR.hh"
#include "ControlFunctions.hh"
#include "WebXRMotorParams.hh"

#define INTERVAL_IN_SEC 10

using namespace ARCS;

ControlFunctions controlFunctions;
long currentTimeInMsec;
double currentTime;

MotorModel MotorModels[ConstParams::ACTUATOR_NUM] = {
    MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
               DEFAULT_VISCOUS_FRICTION_COEFFICIENT),
    MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
               DEFAULT_VISCOUS_FRICTION_COEFFICIENT),
    MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
               DEFAULT_VISCOUS_FRICTION_COEFFICIENT),
    MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
               DEFAULT_VISCOUS_FRICTION_COEFFICIENT),
    MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
               DEFAULT_VISCOUS_FRICTION_COEFFICIENT),
};

extern "C" void invokeFunctionsCallback() {
  currentTimeInMsec += INTERVAL_IN_SEC;
  currentTime = currentTimeInMsec / 1000.0;

  controlFunctions.ControlFunction1(currentTime, 0, 0);
  controlFunctions.ControlFunction2(currentTime, 0, 0);
  controlFunctions.ControlFunction3(currentTime, 0, 0);
  controlFunctions.UpdateControlValue();
  controlFunctions.UpdateScreen();
}

extern "C" void updateModeCallback(int controlFunctionsMode) {
  if (controlFunctionsMode == 1) {
    controlFunctions.UpdateMode(ControlFunctions::CtrlFuncMode::CTRL_INIT);
    currentTimeInMsec = 0;
    currentTime = 0;

    for (int i = 0; i < ConstParams::ACTUATOR_NUM; ++i) {
      MotorModels[i] =
          MotorModel(DEFAULT_TORQUE_CONSTANT, DEFAULT_MOMENT_OF_INERTIA,
                     DEFAULT_VISCOUS_FRICTION_COEFFICIENT);
      MotorModels[i].setPreviousTime(currentTime);
    }
  } else if (controlFunctionsMode == 2) {
    controlFunctions.UpdateMode(ControlFunctions::CtrlFuncMode::CTRL_LOOP);
  } else if (controlFunctionsMode == 3) {
    controlFunctions.UpdateMode(ControlFunctions::CtrlFuncMode::CTRL_EXIT);
  }
}
