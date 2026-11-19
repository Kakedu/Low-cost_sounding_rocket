# Extended Kalman filters used by sil.py.
# FlightKF: state [h, v, Cd].  FlightKF4: adds the static-port coefficient.

import math

import numpy as np

G0 = 9.80665


from atmosphere_model import troposphere_density_kg_m3

def isa_rho(h_m, dT=0.0):
    return troposphere_density_kg_m3(h_m, dT)


def isa_drho_dh(h_m, dh=1.0, dT=0.0):
    return (isa_rho(h_m + dh, dT) - isa_rho(h_m - dh, dT)) / (2 * dh)


class FlightKF:
    def __init__(self, dt, sensors, cd0, dcd_brake, s_ref, mass_kg):
        self.dt = dt
        self.S, self.m, self.dcd = s_ref, mass_kg, dcd_brake
        self.r_baro = sensors["baro_sigma_m"] ** 2
        self.port_cp_sd = sensors.get("port_cp_sd", 0.0)
        self.dT = sensors.get("dT_K", 0.0)
        self.r_acc = (sensors["acc_sigma_g"] * G0) ** 2
        self.x = np.array([0.0, 0.0, cd0])
        self.P = np.diag([1.0, 1.0, 0.1 ** 2])
        self.coasting = False
        self.decel_count = 0
        self.elapsed = 0.0
        self.coast_min_time_s = sensors.get("coast_min_time_s", 0.0)
        self.coast_detected_time_s = None
        self.innovation_nis_limit = sensors.get("innovation_nis_limit", 9.0)
        self.brake_relative_sd = sensors.get("brake_relative_sd", 0.20)
        self.measurement_updates = {"baro": 0, "accel": 0}
        self.measurement_rejections = {"baro": 0, "accel": 0}
        self.max_nis = {"baro": 0.0, "accel": 0.0}
        self.last_innovation = {}
        self.h, self.v, self.cd = 0.0, 0.0, None
        self.cd_projection_count = 0

    def _project_drag(self):
        """Project the mean to Cd>=0; P remains an unconstrained approximation.

        This is a projected EKF, not a truncated-posterior Bayesian update.
        The controller's separate [0.2,1] predictor clamp is unchanged.
        """
        if self.x[2] < 0.0:
            self.x[2] = 0.0
            self.cd_projection_count += 1

    def step(self, a_meas, h_meas, truth, brake):
        # This clock is advanced from the configured sensor period. No plant
        # thrust, propellant contents or true burnout is used by the estimator.
        self.elapsed += self.dt
        phase_ready = self.elapsed + 1e-10 >= self.coast_min_time_s
        self.decel_count = self.decel_count + 1 if phase_ready and a_meas < -0.2 * G0 else 0
        if not self.coasting and self.decel_count >= 3:
            self.coasting = True
            self.coast_detected_time_s = self.elapsed
        if not self.coasting:
            self._predict_boost(a_meas)
            self._update_baro(h_meas)
        else:
            self._predict_coast(brake)
            self._update_baro(h_meas)
            self._update_accel(a_meas, brake)
        self.h, self.v = float(self.x[0]), float(self.x[1])
        self.cd = float(self.x[2]) if self.coasting else None

    def _predict_boost(self, a_meas):
        dt = self.dt
        self.x[0] += self.x[1] * dt
        self.x[1] += (a_meas - G0) * dt
        F = np.array([[1, dt, 0], [0, 1, 0], [0, 0, 1]])
        q = (0.5 * G0) ** 2 * dt
        Q = np.diag([0.0, q, 0.0])
        self.P = F @ self.P @ F.T + Q

    def _drag_accel(self, x, brake):
        h, v, cd = x
        return -isa_rho(h, self.dT) * v * abs(v) * self.S * (cd + self.dcd * brake) / (2 * self.m)

    def _predict_coast(self, brake):
        self._project_drag()
        dt = self.dt
        h, v, cd = self.x
        ad = self._drag_accel(self.x, brake)
        rho = isa_rho(h, self.dT)
        cdt = cd + self.dcd * brake
        dadh = -isa_drho_dh(h, dT=self.dT) * v * abs(v) * self.S * cdt / (2 * self.m)
        dadv = -rho * 2 * abs(v) * self.S * cdt / (2 * self.m)
        dadc = -rho * v * abs(v) * self.S / (2 * self.m)
        self.x = np.array([h + v * dt, v + (ad - G0) * dt, cd])
        F = np.array([[1, dt, 0], [dadh * dt, 1 + dadv * dt, dadc * dt], [0, 0, 1]])
        Q = np.diag([0.0, (0.05) ** 2 * dt, (0.01) ** 2 * dt])
        self.P = F @ self.P @ F.T + Q

    def _update_baro(self, h_meas):
        H = np.array([[1.0, 0.0, 0.0]])
        y = h_meas - self.x[0]
        r = self.r_baro
        if self.port_cp_sd:
            r = r + (self.port_cp_sd * self.x[1] ** 2 / (2 * G0)) ** 2
        self._measurement_update("baro", H, y, r)

    def _measurement_update(self, label, H, innovation, variance):
        """Acceleration-only scalar 3-sigma test, followed by Joseph.

        Rejecting an incompatible acceleration preserves propagated x/P.
        Barometer NIS is recorded without hard rejection: direct paired
        diagnostics found that barometer rejection prevented reacquisition
        and locked a wrong h/v/Cp branch. No h/v clipping or flight dropping.
        """
        S = float((H @ self.P @ H.T).item() + variance)
        if not math.isfinite(S) or S <= 0.0:
            raise ValueError("Invalid innovation covariance")
        nis = float(innovation * innovation / S)
        self.measurement_updates[label] += 1
        self.max_nis[label] = max(self.max_nis[label], nis)
        accepted = label == "baro" or nis <= self.innovation_nis_limit
        self.last_innovation[label] = dict(innovation=float(innovation), S=S,
                                           R=float(variance), NIS=nis, accepted=accepted)
        if not accepted:
            self.measurement_rejections[label] += 1
            return
        Kg = self.P @ H.T / S
        self.x = self.x + (Kg * innovation).ravel()
        A = np.eye(len(self.x)) - Kg @ H
        self.P = A @ self.P @ A.T + variance * Kg @ Kg.T
        self.P = 0.5 * (self.P + self.P.T)
        self._project_drag()

    def _accel_variance(self, brake):
        h, v = self.x[:2]
        brake_specific_force = -isa_rho(h, self.dT) * v * abs(v) * self.S * self.dcd * brake / (2 * self.m)
        # The assigned brake coefficient has relative SD .20. Treating its
        # per-flight fixed error as local model variance is an approximation;
        # it does not represent its exact cross-time covariance.
        return (self.r_acc + (0.02 * G0) ** 2
                + (self.brake_relative_sd * brake_specific_force) ** 2)

    def _update_accel(self, a_meas, brake):
        h, v, cd = self.x
        cdt = cd + self.dcd * brake
        rho = isa_rho(h, self.dT)
        H = np.array([[-isa_drho_dh(h, dT=self.dT) * v * abs(v) * self.S * cdt / (2 * self.m),
                       -rho * 2 * abs(v) * self.S * cdt / (2 * self.m),
                       -rho * v * abs(v) * self.S / (2 * self.m)]])
        y = a_meas - self._drag_accel(self.x, brake)
        self._measurement_update("accel", H, y, self._accel_variance(brake))


# Baro is modelled as h + Cp*v^2/(2g). Cp can be estimated because the
# accelerometer gives v independently.
class FlightKF4(FlightKF):

    def __init__(self, dt, sensors, cd0, dcd_brake, s_ref, mass_kg):
        super().__init__(dt, dict(sensors, port_cp_sd=0.0), cd0, dcd_brake, s_ref, mass_kg)
        sd = sensors.get("port_cp_sd", 0.03) or 0.03
        self.x = np.array([0.0, 0.0, cd0, 0.0])
        self.P = np.diag([1.0, 1.0, 0.1 ** 2, sd ** 2])
        self.cp = 0.0

    def step(self, a_meas, h_meas, truth, brake):
        super().step(a_meas, h_meas, truth, brake)
        self.cp = float(self.x[3])

    def _predict_boost(self, a_meas):
        dt = self.dt
        self.x[0] += self.x[1] * dt
        self.x[1] += (a_meas - G0) * dt
        F = np.eye(4); F[0, 1] = dt
        Q = np.diag([0.0, (0.5 * G0) ** 2 * dt, 0.0, 0.0])
        self.P = F @ self.P @ F.T + Q

    def _predict_coast(self, brake):
        self._project_drag()
        dt = self.dt
        h, v, cd, cp = self.x
        ad = self._drag_accel(self.x[:3], brake)
        rho = isa_rho(h, self.dT)
        cdt = cd + self.dcd * brake
        dadh = -isa_drho_dh(h, dT=self.dT) * v * abs(v) * self.S * cdt / (2 * self.m)
        dadv = -rho * 2 * abs(v) * self.S * cdt / (2 * self.m)
        dadc = -rho * v * abs(v) * self.S / (2 * self.m)
        self.x = np.array([h + v * dt, v + (ad - G0) * dt, cd, cp])
        F = np.eye(4)
        F[0, 1] = dt
        F[1, 0], F[1, 1], F[1, 2] = dadh * dt, 1 + dadv * dt, dadc * dt
        Q = np.diag([0.0, 0.05 ** 2 * dt, 0.01 ** 2 * dt, 0.0])
        self.P = F @ self.P @ F.T + Q

    def _update_baro(self, h_meas):
        h, v, cd, cp = self.x
        H = np.array([[1.0, cp * v / G0, 0.0, v * v / (2 * G0)]])
        y = h_meas - (h + cp * v * v / (2 * G0))
        self._measurement_update("baro", H, y, self.r_baro)

    def _update_accel(self, a_meas, brake):
        h, v, cd, cp = self.x
        cdt = cd + self.dcd * brake
        rho = isa_rho(h, self.dT)
        H = np.array([[-isa_drho_dh(h, dT=self.dT) * v * abs(v) * self.S * cdt / (2 * self.m),
                       -rho * 2 * abs(v) * self.S * cdt / (2 * self.m),
                       -rho * v * abs(v) * self.S / (2 * self.m), 0.0]])
        y = a_meas - self._drag_accel(self.x[:3], brake)
        self._measurement_update("accel", H, y, self._accel_variance(brake))
