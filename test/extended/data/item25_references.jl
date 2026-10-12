# Reference data and conventions are documented in ../README.md.
# ITEM25_ROWS: independent 384-bit Teukolsky integration references near the positive imaginary axis
# (omega = eps + 0.3i). incidence/reflection are B/B_trans (IN) or C_up/C_trans, C_ref/C_trans (UP).
# meas_*: historical package errors at these points, not used as tolerances.
const ITEM25_ROWS = [
    (id = "lower_e10", a = 0.0, omega = complex(1e-10, -0.3), branch = "IN", lambda = complex(4.0, 0.0), incidence = complex(1417.323255170568, -1289.2811272176389), reflection = complex(3.5608563375801725, -7.333711560290893e-10), meas_before = (4.53766334162653e-15, 3.3574624556156806e-15)),
    (id = "lower_e10", a = 0.0, omega = complex(1e-10, -0.3), branch = "UP", lambda = complex(4.0, 0.0), incidence = complex(-531.4962210918634, 483.48042226370103), reflection = complex(2752.606451659051, -2503.6783074265736), meas_before = (3.027564813197553e-14, 4.7123218772514426e-14)),
    (id = "upper_e10", a = 0.0, omega = complex(1e-10, 0.3), branch = "IN", lambda = complex(4.0, 0.0), incidence = complex(3846.3869020421307, 4.5425481196032135e-06), reflection = complex(45.67541400622709, 41.54482106813552), meas_before = (1.0001397386606021, 0.7578128973634095)),
    (id = "upper_e10", a = 0.0, omega = complex(1e-10, 0.3), branch = "UP", lambda = complex(4.0, 0.0), incidence = complex(360.5987720664498, 3.5073914203229084e-07), reflection = complex(214.59326356735198, 7.213446066030415e-07), meas_before = (1.2795748688195132, 0.8792861822595878)),
    (id = "upper_e5", a = 0.0, omega = complex(1e-05, 0.3), branch = "IN", lambda = complex(4.0, 0.0), incidence = complex(3846.386864283189, 0.454254809349808), reflection = complex(45.67192658644627, 41.53999760331278), meas_before = (1.4282977819298787e-06, 2.147880137783631e-10)),
    (id = "upper_e5", a = 0.0, omega = complex(1e-05, 0.3), branch = "UP", lambda = complex(4.0, 0.0), incidence = complex(360.5987695076713, 0.035073914044505425), reflection = complex(214.5932460063031, 0.07213445689717007), meas_before = (1.3235711962439875e-10, 9.096060279614135e-11)),
    (id = "upper_e8", a = 0.0, omega = complex(1e-08, 0.3), branch = "IN", lambda = complex(4.0, 0.0), incidence = complex(3846.386902042093, 0.0004542548119603187), reflection = complex(45.67541055355921, 41.544816292593836), meas_before = (0.3581184082757427, 0.00012078533305802692)),
    (id = "upper_e8", a = 0.0, omega = complex(1e-08, 0.3), branch = "UP", lambda = complex(4.0, 0.0), incidence = complex(360.5987720664472, 3.507391420322893e-05), reflection = complex(214.59326356733442, 7.213446066030039e-05), meas_before = (0.00017955481065194898, 0.0001233847955707206)),
    (id = "upper_spin07", a = 0.7, omega = complex(1e-10, 0.3), branch = "IN", lambda = complex(3.988339448934813, -1.3994548196352707), incidence = complex(2982.3054503576254, 1029.1374158888125), reflection = complex(7.788888043692235, 47.03243752440666), meas_before = (6.36773137193987e-12, 7.392458020053595e-12)),
    (id = "upper_spin07", a = 0.7, omega = complex(1e-10, 0.3), branch = "UP", lambda = complex(3.988339448934813, -1.3994548196352707), incidence = complex(328.9620775709849, -61.79072442323365), reflection = complex(-11.298135410551616, -9.171010618653563), meas_before = (7.39810097524337e-12, 2.4679341769033163e-12)),
]
# ITEM25_AXIS: incidence-type amplitude (B_inc/B_trans for IN, C_up/C_trans for UP) on and near the
# positive imaginary axis (independent 384-bit Teukolsky integration, Wronskian-based).
const ITEM25_AXIS = [
    (a = 0.0, omega = complex(0.0, 0.3), branch = "IN", incidence = complex(3846.3869020421307, -0.0), meas_before = 4.074087471532795e-11),
    (a = 0.0, omega = complex(0.0, 0.3), branch = "UP", incidence = complex(360.5987720664498, -0.0), meas_before = 4.5121916049362776e-15),
    (a = 0.0, omega = complex(1e-10, 0.3), branch = "IN", incidence = complex(3846.3869020421307, 4.5425481196032135e-06), meas_before = 1.0001397386606021),
    (a = 0.0, omega = complex(1e-10, 0.3), branch = "UP", incidence = complex(360.5987720664498, 3.5073914203229084e-07), meas_before = 1.2795748688195132),
    (a = 0.0, omega = complex(1e-08, 0.3), branch = "IN", incidence = complex(3846.386902042093, 0.0004542548119603187), meas_before = 0.3581184082757427),
    (a = 0.0, omega = complex(1e-08, 0.3), branch = "UP", incidence = complex(360.5987720664472, 3.507391420322893e-05), meas_before = 0.00017955481065194898),
    (a = 0.7, omega = complex(0.0, 0.3), branch = "IN", incidence = complex(2982.3054518676413, 1029.1374122806278), meas_before = 1.6327358836554718e-11),
    (a = 0.7, omega = complex(0.0, 0.3), branch = "UP", incidence = complex(328.9620775149858, -61.790724756436724), meas_before = 1.4703400744233603e-11),
    (a = 0.7, omega = complex(1e-10, 0.3), branch = "IN", incidence = complex(2982.3054503576254, 1029.1374158888125), meas_before = 6.36773137193987e-12),
    (a = 0.7, omega = complex(1e-10, 0.3), branch = "UP", incidence = complex(328.9620775709849, -61.79072442323365), meas_before = 7.39810097524337e-12),
    (a = 0.7, omega = complex(1e-08, 0.3), branch = "IN", incidence = complex(2982.305300866095, 1029.137773099064), meas_before = 1.776119126995466e-11),
    (a = 0.7, omega = complex(1e-08, 0.3), branch = "UP", incidence = complex(328.96208311489545, -61.79069143613488), meas_before = 1.7174128132776872e-11),
]

