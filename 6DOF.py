import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ==========================================
# 1. QUATERNION ALGEBRA & KINEMATICS
# ==========================================

def quat_normalize(q):
    norm = np.linalg.norm(q)
    return q / norm if norm > 0 else np.array([1.0, 0.0, 0.0, 0.0])

def quat_conjugate(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def quat_multiply(q1, q2):
    """Hamilton product: q1 * q2"""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def quat_to_dcm(q):
    """Direction Cosine Matrix: Transforms Inertial -> Body (B_C_I)"""
    qw, qx, qy, qz = q
    return np.array([
        [1 - 2*(qy**2 + qz**2),     2*(qx*qy + qw*qz),     2*(qx*qz - qw*qy)],
        [    2*(qx*qy - qw*qz), 1 - 2*(qx**2 + qz**2),     2*(qy*qz + qw*qx)],
        [    2*(qx*qz + qw*qy),     2*(qy*qz - qw*qx), 1 - 2*(qx**2 + qy**2)]
    ])

def quat_derivative(q, omega_b):
    """Kinematic rate equation: dq/dt = 0.5 * Omega(omega_b) * q"""
    p, q_y, r = omega_b
    Omega = np.array([
        [ 0,  -p, -q_y,  -r],
        [ p,   0,    r, -q_y],
        [ q_y, -r,   0,   p],
        [ r,  q_y,  -p,   0]
    ])
    return 0.5 * Omega @ q

# ==========================================
# 2. ROCKET CONFIGURATION & CONTROLLER GAINS
# ==========================================

class TVCRocket:
    def __init__(self):
        # Mass & Inertia
        self.mass = 50.0  # kg
        self.I_body = np.diag([0.2, 35.0, 35.0])
        self.I_inv = np.linalg.inv(self.I_body)
        
        # Geometry & Aerodynamics
        self.reference_area = np.pi * (0.15 / 2.0)**2
        self.cop_offset_x = 0.25   # Distance of CP behind CG (m)
        self.Cd_axial = 0.40
        self.C_normal_alpha = 6.0
        self.damping_pitch = 10.0
        
        # Propulsion & TVC Gimbal
        self.thrust_magnitude = 2400.0  # N (~2.4 kN)
        self.burn_time = 4.0            # s
        self.gimbal_arm = 1.2           # Gimbal hinge distance behind CG (m)
        self.max_gimbal_angle = np.deg2rad(5.0)  # Max deflection ±5 deg
        
        # PD Attitude Controller Gains (Pitch & Yaw)
        self.Kp = 85.0   # Proportional gain
        self.Kd = 22.0   # Derivative gain (damping)

def atmospheric_density(alt):
    if alt < 0: alt = 0
    return 1.225 * np.exp(-alt / 8500.0)

# ==========================================
# 3. TVC CONTROLLER LOGIC
# ==========================================

def compute_tvc_deflection(q_current, omega_b, q_target, rocket):
    """
    Computes gimbal angles [delta_y, delta_z] via quaternion error feedback.
    q_err = q_target^-1 * q_current
    """
    q_err = quat_multiply(quat_conjugate(q_target), q_current)
    
    # Enforce shortest rotation path
    if q_err[0] < 0:
        q_err = -q_err
        
    # Attitude error vector (small angle approx: 2 * vector part)
    # err_b[1] -> pitch error, err_b[2] -> yaw error
    err_b = 2.0 * q_err[1:4]
    
    # PD feedback torque commands in body frame
    # M = -Kp * e - Kd * omega
    tau_cmd_pitch = -(rocket.Kp * err_b[2] + rocket.Kd * omega_b[1])
    tau_cmd_yaw   = -(rocket.Kp * err_b[1] + rocket.Kd * omega_b[2])
    
    # TVC torque relation:
    # Gimbal deflection delta_pitch (around body Y) creates a force in body Z:
    # tau_pitch = -gimbal_arm * F_thrust * sin(delta_pitch)
    # For small angles: delta ~ -tau / (arm * Thrust)
    denom = rocket.gimbal_arm * rocket.thrust_magnitude
    
    delta_pitch = np.clip(tau_cmd_pitch / denom, -rocket.max_gimbal_angle, rocket.max_gimbal_angle)
    delta_yaw   = np.clip(-tau_cmd_yaw / denom,  -rocket.max_gimbal_angle, rocket.max_gimbal_angle)
    
    return delta_pitch, delta_yaw

# ==========================================
# 4. 6DOF DYNAMICS FUNCTION
# ==========================================

def rocket_6dof_tvc_dynamics(t, state, rocket, q_target):
    r_i = state[0:3]
    v_i = state[3:6]
    q   = quat_normalize(state[6:10])
    w_b = state[10:13]
    
    alt = r_i[2]
    rho = atmospheric_density(alt)
    
    B_C_I = quat_to_dcm(q)
    I_C_B = B_C_I.T
    v_b = B_C_I @ v_i
    v_mag = np.linalg.norm(v_b)
    
    F_b = np.zeros(3)
    M_b = np.zeros(3)
    
    # --- 1. TVC Propulsion ---
    if t <= rocket.burn_time:
        delta_p, delta_y = compute_tvc_deflection(q, w_b, q_target, rocket)
        
        # Thrust unit vector deflected by gimbal angles
        # delta_p pitches nozzle up/down (Z axis), delta_y turns left/right (Y axis)
        thrust_dir = np.array([
            np.cos(delta_p) * np.cos(delta_y),
            np.sin(delta_y),
            np.sin(delta_p)
        ])
        thrust_force_b = rocket.thrust_magnitude * thrust_dir
        F_b += thrust_force_b
        
        # Gimbal torque: tau = r_gimbal x F_thrust
        r_gimbal = np.array([-rocket.gimbal_arm, 0.0, 0.0])
        M_b += np.cross(r_gimbal, thrust_force_b)
    
    # --- 2. Aerodynamics ---
    if v_mag > 1.0:
        q_dyn = 0.5 * rho * (v_mag**2)
        
        # Axial drag
        F_b[0] -= q_dyn * rocket.reference_area * rocket.Cd_axial
        
        # Lift / Normal aerodynamic restoring forces
        alpha = np.arctan2(-v_b[2], v_b[0])
        beta  = np.arctan2(v_b[1], v_b[0])
        
        F_aero_y = -q_dyn * rocket.reference_area * rocket.C_normal_alpha * beta
        F_aero_z =  q_dyn * rocket.reference_area * rocket.C_normal_alpha * alpha
        F_b[1] += F_aero_y
        F_b[2] += F_aero_z
        
        # CP moments & aero damping
        M_b[1] += rocket.cop_offset_x * F_aero_z - rocket.damping_pitch * w_b[1]
        M_b[2] -= rocket.cop_offset_x * F_aero_y + rocket.damping_pitch * w_b[2]
        M_b[0] -= 0.05 * w_b[0]
        
    # --- 3. Translation & Rotation Integration ---
    g_i = np.array([0.0, 0.0, -9.80665])
    a_i = (1.0 / rocket.mass) * (I_C_B @ F_b) + g_i
    
    dq_dt = quat_derivative(q, w_b)
    dw_b_dt = rocket.I_inv @ (M_b - np.cross(w_b, rocket.I_body @ w_b))
    
    return np.hstack([v_i, a_i, dq_dt, dw_b_dt])

# ==========================================
# 5. SIMULATION & TRAJECTORY PLOT
# ==========================================

def ground_event(t, state, *args):
    return state[2]
ground_event.terminal = True
ground_event.direction = -1

if __name__ == "__main__":
    rocket = TVCRocket()
    
    # Command a gravity turn / pitch program:
    # Target: 80 degrees pitch (10 degrees off vertical along +X)
    target_pitch_deg = 80.0
    theta_target = np.deg2rad(90.0 - target_pitch_deg)
    q_target = np.array([np.cos(theta_target / 2.0), 0.0, np.sin(theta_target / 2.0), 0.0])
    
    # Initial state: Purely vertical on launch rail (pitch = 90 deg) with a slight disturbance
    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    r0 = np.array([0.0, 0.0, 0.5])
    v0 = np.array([0.0, 0.0, 0.0])
    w0 = np.array([0.0, 0.05, -0.03])  # Initial tip-off rates (rad/s)
    
    state0 = np.hstack([r0, v0, q0, w0])
    t_span = (0.0, 35.0)
    
    sol = solve_ivp(
        fun=rocket_6dof_tvc_dynamics,
        t_span=t_span,
        y0=state0,
        args=(rocket, q_target),
        method='RK45',
        t_eval=np.linspace(t_span[0], t_span[1], 800),
        events=ground_event,
        rtol=1e-6,
        atol=1e-8
    )
    
    # Output metrics
    apogee_idx = np.argmax(sol.y[2])
    print(f"Apogee reached: {sol.y[2, apogee_idx]:.2f} m at t = {sol.t[apogee_idx]:.2f} s")
    print(f"Max velocity:   {np.max(np.linalg.norm(sol.y[3:6], axis=0)):.2f} m/s")
    print(f"Downrange dist: {sol.y[0, -1]:.2f} m")