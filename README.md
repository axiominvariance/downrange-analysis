# Rocket Flight Safety Analysis

## Monte Carlo Simulation with Two-Stage Rocket Dynamics

A physics-based Monte Carlo simulation for analyzing rocket launch safety and debris hazard zones. Developed independently, outside of university coursework, as part of a self-directed portfolio in computational physics.

---

## Project structure

```
monte-carlo-rocket/
├── main.py                 Entry point: loads config, runs simulation, plots results
├── rocket_simulator.py     Core physics engine (StageParameters, RocketSimulator)
├── config.json             All simulation parameters (editable)
├── test_physics.py         Unit tests
├── Herleitung.pdf          Handwritten derivation of the equations of motion (German)
├── requirements.txt        Python dependencies
└── README.md
```

`main.py` handles configuration loading, visualization (matplotlib), and the top-level simulation flow. `rocket_simulator.py` contains the physics: the equations of motion, the ODE integration (via `scipy.integrate.solve_ivp`), stage separation logic, and the Monte Carlo sampling loop. `config.json` stores all physical parameters and uncertainty definitions so the simulation can be reconfigured without touching code.

---

## Physics model

### State vector

The rocket's motion is described by four coupled first-order ODEs. The state vector is:

```
y = [x, vx, y, vy]
```

where `x` and `y` are horizontal and vertical position (m), and `vx`, `vy` are the corresponding velocities (m/s).

### Forces

Three forces act on the rocket at every instant:

**Thrust:**

$$\vec{F}_{\text{thrust}} = F \cdot (\cos\theta,\; \sin\theta)$$

During Stage 1, the thrust angle θ is fixed at the launch angle (85°). The thrust magnitude F is constant during a burn and zero after burnout.

**Gravity:**

$$\vec{F}_g = (0,\; -m(t) \cdot g)$$

**Drag (quadratic model):**

$$\vec{F}_{\text{drag}} = -\tfrac{1}{2}\,\rho\, C_d\, A\, |\vec{v}_{\text{rel}}| \cdot \vec{v}_{\text{rel}}, \qquad \vec{v}_{\text{rel}} = \vec{v}_{\text{rocket}} - \vec{v}_{\text{wind}}$$

Wind is modeled as a constant horizontal velocity. The drag force opposes the velocity relative to the air mass.

### Mass evolution (Tsiolkovsky)

During a burn, mass decreases linearly:

$$m(t) = m_0 - \dot{m} \cdot t, \qquad \dot{m} = \frac{m_{\text{propellant}}}{t_{\text{burn}}}$$

After burnout, mass is constant at $m_{\text{dry}} = m_{\text{total}} - m_{\text{propellant}}$.

### Equations of motion

$$\dot{x} = v_x$$

$$\dot{v}_x = \frac{1}{m(t)}\left(F\cos\theta - \tfrac{1}{2}\rho\, C_d\, A\, |v_{\text{rel}}|\, v_{x,\text{rel}}\right)$$

$$\dot{y} = v_y$$

$$\dot{v}_y = \frac{1}{m(t)}\left(F\sin\theta - m(t)\,g - \tfrac{1}{2}\rho\, C_d\, A\, |v_{\text{rel}}|\, v_{y,\text{rel}}\right)$$

These are integrated numerically using `scipy.integrate.solve_ivp` with the RK45 method.

### Stage separation (multi-body)

At $t = t_{b1}$ (end of Stage 1 burn), the state vector $(x, v_x, y, v_y)$ at separation serves as initial conditions for two independent trajectories:

1. **Stage 2 + payload**: Continues powered flight with Stage 2 thrust parameters, then transitions to ballistic flight after Stage 2 burnout.
2. **Stage 1 debris**: Empty first stage ($m = m_{\text{dry},1}$, thrust = 0) follows a ballistic arc until ground impact.

### Ballistic flight

After burnout, thrust is zero and mass is constant. The same ODEs apply with $F_{\text{thrust}} = 0$.

---

## Monte Carlo method

Each simulation parameter has an associated uncertainty modeled as a Gaussian:

| Parameter | Nominal | 1σ uncertainty |
|---|---|---|
| Stage 1 thrust | 15,000 N | ±3% |
| Stage 2 thrust | 5,000 N | ±3% |
| Stage 1 total mass | 800 kg | ±1% |
| Wind speed | 0 m/s | ±5 m/s |
| Launch angle | 85° | ±0.5° |

The simulation runs 10,000 times with randomly sampled parameters. The 95% confidence interval is the 2.5th and 97.5th percentile of the resulting impact distributions.

---

## Output

### Nominal trajectory

Payload impact: ~5.51 km downrange. Stage 1 debris impact: ~2.99 km downrange. Maximum altitude: ~37 km.

### Monte Carlo results (10,000 samples)

```
Payload impact:   mean 5.5 km,  95% CI: 1.5 – 9.5 km
Stage 1 debris:   mean 3.0 km,  95% CI: 1.4 – 4.6 km
```

---

## Known shortcomings

### Constant air density

This is the most significant physics simplification. The simulation uses a constant sea-level air density (ρ = 1.225 kg/m³) at all altitudes. In reality, air density decreases roughly exponentially with altitude — at 37 km (the peak altitude of this rocket), the actual density is about 0.5% of the sea-level value. This means drag forces are overestimated by a factor of ~200 at peak altitude, which compresses the trajectory and drastically underestimates the true downrange impact distance. The correct approach is to use the barometric height formula (barometrische Höhenformel), where ρ(h) follows the International Standard Atmosphere: a linear temperature lapse in the troposphere (0–11 km) transitioning to an isothermal model in the lower stratosphere (11–25 km), with density computed from the ideal gas law at each altitude.

### Fixed thrust direction

The thrust vector angle θ is fixed at the launch angle (85°) for the entire flight — both stages. Real rockets perform a gravity turn during upper-stage flight, where the thrust vector tracks the velocity direction: θ(t) = atan2(vy, vx). This naturally pitches the rocket from near-vertical to near-horizontal, which is the most efficient way to build orbital or downrange velocity. With a fixed angle, the trajectory stays much more vertical than it should, further underestimating downrange distance.

### Mass sampling in Monte Carlo

The Monte Carlo samples `m_total` directly from a Gaussian while keeping `m_propellant` fixed. If the sampled `m_total` falls below `m_propellant` (rare at 1% standard deviation, but possible), the simulation would produce an unphysical negative dry mass. The current code catches this via a bare `except` and marks it as NaN, which hides the problem. The correct approach is to sample `m_dry` (the structurally uncertain quantity) and reconstruct `m_total = m_dry + m_propellant`. This eliminates the edge case entirely and is physically more meaningful — structural mass is what actually varies in manufacturing, not propellant load.

### No convergence check

The Monte Carlo uses a fixed 10,000 samples without verifying that the statistics have converged. A proper analysis would run the simulation at increasing sample sizes (e.g. 1k, 5k, 10k, 20k) and show that the 95% CI stabilizes. Without this, there is no evidence that 10,000 samples is sufficient.

### Independent parameters only

All uncertain parameters are sampled independently. In reality, thrust and specific impulse are correlated (they both depend on combustion chamber conditions), and structural mass components may co-vary. A correlation matrix or copula-based sampling would capture these dependencies. The current model likely overestimates the spread of the impact distribution because it allows physically unlikely parameter combinations.

### Python performance

The Monte Carlo loop runs 10,000 full trajectory integrations sequentially in Python. Each integration calls `scipy.integrate.solve_ivp`, which is implemented in compiled code, but the outer loop, parameter sampling, and object construction are pure Python. This makes the simulation slow — on the order of minutes for 10,000 samples. A C++ or Rust implementation with a hand-written RK4 integrator would run the same Monte Carlo in 1–2 seconds. For production use or larger sample sizes, Python is not the right language for this kind of compute-bound simulation. Additionally, the Monte Carlo loop is embarrassingly parallel (each run is independent), but the current implementation is single-threaded.

### No validation against analytical solutions

The simulation has unit tests, but no comparison against a known analytical solution. Even a simple check — vertical launch with no drag should give $h_{\max} = v^2 / (2g)$ after burnout — would provide confidence that the integrator and force model are correct. Without this, the simulation is only tested for internal consistency, not physical correctness.

---

## Dependencies

- Python 3.x
- NumPy ≥ 1.20
- SciPy ≥ 1.7
- Matplotlib ≥ 3.4

---

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python main.py
```

---

## Author

Silas-Anton Voglmaier — B.Sc. Physics, University of Augsburg
GitHub: [@axiominvariance](https://github.com/axiominvariance)