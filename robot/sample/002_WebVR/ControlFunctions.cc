#include "ControlFunctions.hh"
#include "ARCSscrparams.hh"
#include "ArcsMatrix.hh"
#include <math.h>

using namespace ARCS;
using namespace ArcsMatrix;

namespace {
int currentStep = 1;

ArcsMat<ConstParams::ACTUATOR_NUM, 1> PositionRes;
ArcsMat<ConstParams::ACTUATOR_NUM, 1> CurrentRef;
} // namespace

bool ControlFunctions::ControlFunction1(double t, double Tact, double Tcmp) {
  [[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[0] * 1e-9;

  if (CmdFlag == CTRL_INIT) {
    Initializing = true;
    Screen.InitOnlineSetVar();
    Interface.ServoON();
    Initializing = false;

    for (unsigned int i = 0; i < ConstParams::ACTUATOR_NUM; i += 1) {
      PositionRes(i, 1) = 0;
    }
  }

  if (CmdFlag == CTRL_LOOP) {
    Interface.GetPosition(PositionRes);

    double OnlineSetVar = 0.0;
    Screen.GetOnlineSetVar(OnlineSetVar);
    Screen.SetVarIndicator(OnlineSetVar);

    Graph.SetTime(Tact, t);

    for (int i = 0; i < 4; i += 1) {
      Graph.SetVars(i, PositionRes[i] / 2.0 / M_PI, CurrentRef[i] / 2.0 / M_PI);
    }

    if (currentStep == 1) {
      int axisNumber = 1;
      double commandRadian = -M_PI / 2.0;

      CurrentRef = {-3, 0, 0, 0, 0};
      Interface.SetCurrent(CurrentRef);

      double currentRadian = PositionRes(axisNumber, 1);

      if (currentRadian < commandRadian / 1.5) {
        currentStep += 1;
      }
    } else if (currentStep == 2) {
      int axisNumber = 2;
      double commandRadian = -M_PI / 4.0;

      CurrentRef = {0, -3, 0, 0, 0};
      Interface.SetCurrent(CurrentRef);

      double currentRadian = PositionRes(axisNumber, 1);

      if (currentRadian < commandRadian / 1.5) {
        currentStep += 1;
      }
    } else if (currentStep == 3) {
      int axisNumber = 3;
      double commandRadian = (-3 / 4.0) * M_PI;

      CurrentRef = {0, 0, -3, 0, 0};
      Interface.SetCurrent(CurrentRef);

      double currentRadian = PositionRes(axisNumber, 1);

      if (currentRadian < commandRadian / 1.5) {
        currentStep += 1;
      }
    } else if (currentStep == 4) {
      int axisNumber = 1;
      double commandRadian = (1 / 2.0) * M_PI;

      CurrentRef = {3, 0, 0, 0, 0};
      Interface.SetCurrent(CurrentRef);

      double currentRadian = PositionRes(axisNumber, 1);

      if (currentRadian > commandRadian / 1.5) {
        currentStep += 1;
      }
    } else if (currentStep == 5) {
      int axisNumber = 1;
      double commandRadian = (-1 / 2.0) * M_PI;

      CurrentRef = {-3, 0, 0, 0, 0};
      Interface.SetCurrent(CurrentRef);

      double currentRadian = PositionRes(axisNumber, 1);

      if (currentRadian < commandRadian / 1.5) {
        currentStep -= 1;
      }
    }
  }

  if (CmdFlag == CTRL_EXIT) {
    Interface.SetZeroCurrent();
    Interface.ServoOFF();
  }

  return true;
}

bool ControlFunctions::ControlFunction2(double t, double Tact, double Tcmp) {
  [[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[1] * 1e-9;

  if (CmdFlag == CTRL_INIT) {
  }
  if (CmdFlag == CTRL_LOOP) {
  }
  if (CmdFlag == CTRL_EXIT) {
  }

  return true;
}

bool ControlFunctions::ControlFunction3(double t, double Tact, double Tcmp) {
  [[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[2] * 1e-9;

  if (CmdFlag == CTRL_INIT) {
  }
  if (CmdFlag == CTRL_LOOP) {
  }
  if (CmdFlag == CTRL_EXIT) {
  }

  return true;
}

void ControlFunctions::UpdateControlValue() {
  Screen.SetNetworkLink(NetworkLink);
  Screen.SetInitializing(Initializing);
  Screen.SetCurrentAndPosition(CurrentRef, PositionRes);
}
