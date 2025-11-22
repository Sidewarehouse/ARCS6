#ifndef MOTOR_MODEL
#define MOTOR_MODEL

extern "C" void debug(double value);

class MotorModel {
private:
  /** トルク定数 Kt、単位は [Nm/A] です。*/
  double _torqueConstant;

  /** 慣性モーメント J、単位は [kgm²] です。*/
  double _momentOfInertia;

  /** 粘性摩擦係数 D、単位は [Nm/(rad/s)] です。*/
  double _viscousFrictionCoefficient;

  /** 直前のトルク τ、単位は [Nm] です。*/
  double _lastTorque;

  /** 直前の角加速度 α、単位は [rad/s²] です。*/
  double _lastAngularAcceleration;

  /** 直前の角速度 ω、単位は [rad/s] です。*/
  double _lastAngularVelocity;

  /** 直前の位置 θ、単位は [rad] です。*/
  double _lastPosition;

  /** 前回の時刻、単位は [s] です。*/
  double _previousTime;

public:
  /**
   * MotorModel のコンストラクターです。
   *
   * @param torqueConstant トルク定数 Kt、単位は [Nm/A] です。
   * @param momentOfInertia 慣性モーメント J、単位は [kgm²] です。
   * @param viscousFrictionCoefficient 粘性摩擦係数 D、単位は [Nm/(rad/s)]
   * です。
   */
  MotorModel(double torqueConstant, double momentOfInertia,
             double viscousFrictionCoefficient)
      : _torqueConstant(torqueConstant), _momentOfInertia(momentOfInertia),
        _viscousFrictionCoefficient(viscousFrictionCoefficient), _lastTorque(0),
        _lastAngularAcceleration(0), _lastAngularVelocity(0), _lastPosition(0),
        _previousTime(0) {}

  /**
   * 現在の時刻と入力電流から位置を計算します。
   *
   * @param currentTime 現在の時刻、単位は [s] です。
   * @param inputCurrent 入力電流 i、単位は [A] です。
   * @param savePreviousTime true の場合、現在の時刻を保存します。
   * @return 位置 θ、単位は [rad] です。
   */
  double calculatePosition(double currentTime, double inputCurrent,
                           bool savePreviousTime = true) {
    double previousAngularVelocity = _lastAngularVelocity;
    double currentAngularVelocity =
        calculateAngularVelocity(currentTime, inputCurrent);
    double positionDelta =
        calculateDefiniteIntegralValue(_previousTime, previousAngularVelocity,
                                       currentTime, currentAngularVelocity);

    if (savePreviousTime) {
      setPreviousTime(currentTime);
    }

    _lastPosition += positionDelta;
    return _lastPosition;
  }

  /**
   * 現在の時刻と入力電流から角速度 ω [rad/s] を計算します。
   *
   * @param currentTime 現在の時刻、単位は [s] です。
   * @param inputCurrent 入力電流、単位は [A] です。
   * @return 角速度 ω、単位は [rad/s] です。
   */
  double calculateAngularVelocity(double currentTime, double inputCurrent) {
    double previousAngularAcceleration = _lastAngularAcceleration;
    double currentAngularAcceleration =
        calculateAngularAcceleration(currentTime, inputCurrent);
    double angularVelocityDelta = calculateDefiniteIntegralValue(
        _previousTime, previousAngularAcceleration, currentTime,
        currentAngularAcceleration);

    _lastAngularVelocity += angularVelocityDelta;
    return _lastAngularVelocity;
  }

  /**
   * 現在の時刻と入力電流から角加速度 α [rad/s²] を計算します。
   *
   * @param currentTime 現在の時刻、単位は [s] です。
   * @param inputCurrent 入力電流、単位は [A] です。
   * @return 角加速度 α、単位は [rad/s²] です。
   */
  double calculateAngularAcceleration(double currentTime, double inputCurrent) {
    double torque = calculateTorque(currentTime, inputCurrent);
    double frictionalTorque =
        _viscousFrictionCoefficient * _lastAngularVelocity;
    double loadTorque = getLoadTorque(currentTime);
    double sumOfTorque = torque - (frictionalTorque + loadTorque);

    _lastAngularAcceleration = sumOfTorque / _momentOfInertia;
    return _lastAngularAcceleration;
  }

  /**
   * 入力電流からトルク τ [Nm] を計算します。
   *
   * @param currentTime 現在の時刻、単位は [s] です。
   * @param inputCurrent 入力電流、単位は [A] です。
   * @return トルク τ、単位は [Nm] です。
   */
  double calculateTorque(double currentTime, double inputCurrent) {
    _lastTorque = _torqueConstant * inputCurrent;
    return _lastTorque;
  }

  /**
   * 指定時刻における負荷トルクです。
   * このメソッドはオーバーライドして使用されることを想定しています。
   *
   * @param currentTime 現在の時刻、単位は [s] です。
   * @return 角加速度 α、単位は [rad/s²] です。
   */
  double getLoadTorque(double currentTime) { return 0; }

  /**
   * 定積分値を台形公式を使って計算します。
   *
   * @param previousTime 直前の時刻、単位は [s] です。
   * @param previousValue 直前の値、単位は任意です。
   * @param currentTime 現在の時刻、単位は [s] です。
   * @param currentValue 現在の値、単位は任意です。
   * @return 定積分値です。
   */
  double calculateDefiniteIntegralValue(double previousTime,
                                        double previousValue,
                                        double currentTime,
                                        double currentValue) {
    double upperBase = previousValue;
    double lowerBase = currentValue;
    double height = currentTime - previousTime;

    return (upperBase + lowerBase) * height / 2;
  }

  double getTorqueConstant() const { return _torqueConstant; }

  double getMomentOfInertia() const { return _momentOfInertia; }

  double getViscousFrictionCoefficient() const {
    return _viscousFrictionCoefficient;
  }

  double getLastTorque() const { return _lastTorque; }

  double getLastAngularAcceleration() const { return _lastAngularAcceleration; }

  double getLastAngularVelocity() const { return _lastAngularVelocity; }

  double getLastPosition() const { return _lastPosition; }

  double getPreviousTime() const { return _previousTime; }

  void setPreviousTime(double previousTime) { _previousTime = previousTime; }
};

#endif
