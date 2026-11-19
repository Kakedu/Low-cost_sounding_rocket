# 顶点窗口压力选择的只读敏感性

Pressure-selector sensitivity only. Accepted 16.74 K remains the mean of the original apex-window pressure minima after a consistent geometric/geopotential conversion. No flight/MC reruns or model/source changes.

TeleMetrum: fitted GPS apex 29.84940s; minimum 52425Pa at 29.61s (coast).
- minimum_in_apex_window: 52425.0Pa, dT 16.7055K, n=136
- median_in_apex_window: 52511.0Pa, dT 17.6325K, n=136
- minimum_of_coast_state_in_apex_window: 52425.0Pa, dT 16.7055K, n=116
- median_of_coast_state_in_apex_window: 52509.0Pa, dT 17.6109K, n=116
TeleMega: fitted GPS apex 29.52683s; minimum 52393Pa at 28.81s (drogue), 28.82s (drogue).
- minimum_in_apex_window: 52393.0Pa, dT 16.7699K, n=30
- median_in_apex_window: 52417.5Pa, dT 17.0341K, n=30

CSV state_name is the recorder state, not independently measured physical charge/deployment timing.
The TeleMetrum minimum is close to the plotted 29.6 s pressure-disturbance marker; the TeleMega minimum has recorder state drogue. Selection/ejection/static-port behaviour can contribute to the inferred column offset.
A minimum-pressure selector is an apex-window diagnostic, not an independent mean air-temperature observation. Correcting the height convention does not validate its accuracy.
Median and coast-state alternatives are post-hoc sensitivity calculations, not replacement accepted calibration or uncertainty bounds.
