#ifndef INTERFACEFUNCTIONS
#define INTERFACEFUNCTIONS

#include "ArcsMatrix.hh"
#include "ConstParams.hh"
#include "WebXR.hh"
#include <array>

namespace ARCS {
class InterfaceFunctions {
public:
  void ServoON() {}
  void ServoOFF() {}
  void SetZeroCurrent() {}

  void SetCurrent(const ArcsMat<ConstParams::ACTUATOR_NUM, 1> &CurrentRef) {
    for (int i = 1; i <= ConstParams::ACTUATOR_NUM; ++i) {
      double position =
          MotorModels[i - 1].calculatePosition(currentTime, CurrentRef(i, 1));
      setAxisRadian(i, position);
    }
  }

  void GetPosition(ArcsMat<ConstParams::ACTUATOR_NUM, 1> &Position) {
    for (int i = 1; i <= ConstParams::ACTUATOR_NUM; ++i) {
      Position(i, 1) = MotorModels[i - 1].getLastPosition();
    }
  }
};
} // namespace ARCS

#endif
