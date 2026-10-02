# 6DOF Rocket Flight Dynamics Simulator

A 6-Degree-of-Freedom (6DOF) rigid-body rocket trajectory simulator written in Python. It models full 3D spatial translation and rotational dynamics using singularity-free quaternion kinematics, passive aerodynamic restoring stability, and an active Thrust Vector Control (TVC) attitude stabilization system.

---

## Simulated Vehicle Class

The default configuration models a **50 kg High-Power Sounding Rocket** (L/M-class motor equivalent, typical of collegiate competitions such as Spaceport America Cup):
- **Wet Mass:** 50.0 kg
- **Diameter:** 15 cm (0.15 m)
- **Thrust:** 2,400 N (~2.4 kN)
- **Burn Time:** 4.0 seconds
- **Guidance:** Active 2-axis TVC gimbal ($\pm 5^\circ$ max deflection)

All vehicle parameters can be scaled down to hobby model rockets (e.g., 1 kg TVC models) or up to suborbital vehicles by editing `TVCRocket` properties.

---

## Features

- **Singularity-Free Attitude Kinematics:** Uses Hamilton quaternions ($\mathbf{q} = [q_w, q_x, q_y, q_z]$) to propagate orientation without suffering from gimbal lock during vertical flight.
- **Full Rigid-Body Dynamics:** Evaluates Newton-Euler equations of motion, including gyroscopic cross-coupling ($\boldsymbol{\omega} \times \mathbf{I}\boldsymbol{\omega}$) and full inertia tensors.
- **Aerodynamic Modeling:** Models dynamic pressure, exponential atmospheric density decay, axial drag, angle-of-attack normal forces, and Center-of-Pressure (CP) restoring moments.
- **Thrust Vector Control (TVC):** Implements a Proportional-Derivative (PD) quaternion-error feedback controller calculating dynamic nozzle gimbal deflection angles ($\delta_{\text{pitch}}, \delta_{\text{yaw}}$).
- **Adaptive Numerical Integration:** Integrated via SciPy's `solve_ivp` using an adaptive Runge-Kutta (RK45) scheme with automated ground-impact event detection.

---

## Mathematical Formulation

### 1. State Vector (13 States)
$$\mathbf{x}(t) = \big[\, r_x, r_y, r_z, \; v_x, v_y, v_z, \; q_w, q_x, q_y, q_z, \; p, q_y, r \,\big]^T$$

- $\mathbf{r}_i$: Position in Inertial Frame ($+Z$ = Altitude)
- $\mathbf{v}_i$: Velocity in Inertial Frame
- $\mathbf{q}$: Rotation quaternion from Inertial to Body Frame
- $\boldsymbol{\omega}_b = [p, q_y, r]^T$: Angular velocity in Body Frame

### 2. Equations of Motion
- **Translational Acceleration:**
  $$\dot{\mathbf{v}}_i = \frac{1}{m} \mathbf{R}(q)^T \mathbf{F}_b + \mathbf{g}_i$$
- **Quaternion Kinematics:**
  $$\dot{\mathbf{q}} = \frac{1}{2} \boldsymbol{\Omega}(\boldsymbol{\omega}_b) \mathbf{q}$$
- **Euler Rotational Dynamics:**
  $$\dot{\boldsymbol{\omega}}_b = \mathbf{I}^{-1} \Big( \mathbf{M}_b - \boldsymbol{\omega}_b \times (\mathbf{I}\boldsymbol{\omega}_b) \Big)$$

