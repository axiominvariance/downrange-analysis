import numpy as np
from numpy.random import default_rng
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Tuple, Dict, Optional


@dataclass
class StageParameters:
    """Physical parameters for a rocket stage."""
    m_total: float          # Total mass at ignition (kg)
    m_propellant: float     # Propellant mass (kg)
    thrust: float           # Thrust force (N)
    burn_time: float        # Burn duration (s)
    Cd: float               # Drag coefficient
    A: float                # Cross-sectional area (m²)
    
    def validate(self) -> None:
        if self.m_total <= 0:
            raise ValueError(f"m_total must be positive, got {self.m_total}")
        if self.m_propellant < 0 or self.m_propellant > self.m_total:
            raise ValueError(
                f"m_propellant must be in [0, {self.m_total}], got {self.m_propellant}"
            )
        if self.thrust < 0:
            raise ValueError(f"thrust must be non-negative, got {self.thrust}")
        if self.burn_time <= 0:
            raise ValueError(f"burn_time must be positive, got {self.burn_time}")
        if self.Cd <= 0:
            raise ValueError(f"Cd must be positive, got {self.Cd}")
        if self.A <= 0:
            raise ValueError(f"A must be positive, got {self.A}")
    
    @property
    def m_dry(self) -> float:
        """Dry mass (without propellant)."""
        return self.m_total - self.m_propellant


class RocketSimulator:
    """
    Two-stage rocket flight simulator with physics-based equations of motion.
    
    State vector: [x, vx, y, vy]
        x, y     — horizontal/vertical position (m)
        vx, vy   — horizontal/vertical velocity (m/s)
    """
    
    def __init__(
        self,
        stage1: StageParameters,
        stage2: StageParameters,
        launch_angle_rad: float,
        g: float = 9.81,
        rho: float = 1.225,
    ):

        # Validate
        stage1.validate()
        stage2.validate()
        
        if not 0 <= launch_angle_rad <= np.pi/2:
            raise ValueError(
                f"launch_angle_rad must be in [0, π/2], got {launch_angle_rad}"
            )
        
        self.stage1 = stage1
        self.stage2 = stage2
        self.theta = launch_angle_rad
        self.g = g
        self.rho = rho
    
    def _rocket_eom(
        self,
        t: float,
        state: np.ndarray,
        stage: StageParameters,
        t_start: float,
        wind_speed: float,
    ) -> list:
        """
        Equations of motion for powered flight.
        
        Args:
            t: Current time (s)
            state: State vector [x, vx, y, vy]
            stage: Stage parameters
            t_start: Stage ignition time (s)
            wind_speed: Horizontal wind speed (m/s)
        
        Returns:
            State derivatives [dx/dt, dvx/dt, dy/dt, dvy/dt]
        """
        x, vx, y, vy = state
        
        # Mass at current time
        t_burn = t - t_start
        if t_burn < stage.burn_time:
            m_dot = stage.m_propellant / stage.burn_time
            m = stage.m_total - m_dot * t_burn
            F_thrust = stage.thrust
        else:
            m = stage.m_dry
            F_thrust = 0.0
        
        # Relative velocity (rocket - wind)
        vx_rel = vx - wind_speed
        vy_rel = vy
        v_rel_mag = np.sqrt(vx_rel**2 + vy_rel**2)
        
        # Avoid division by zero in drag calculation
        if v_rel_mag < 1e-6:
            Fx_drag = 0.0
            Fy_drag = 0.0
        else:
            Fx_drag = -0.5 * self.rho * stage.Cd * stage.A * v_rel_mag * vx_rel
            Fy_drag = -0.5 * self.rho * stage.Cd * stage.A * v_rel_mag * vy_rel
        
        # Accelerations
        ax = (F_thrust * np.cos(self.theta) + Fx_drag) / m
        ay = (F_thrust * np.sin(self.theta) - m * self.g + Fy_drag) / m
        
        return [vx, ax, vy, ay]
    
    def _ballistic_eom(
        self,
        t: float,
        state: np.ndarray,
        Cd: float,
        A: float,
        m: float,
        wind_speed: float,
    ) -> list:
    
        x, vx, y, vy = state
        
        vx_rel = vx - wind_speed
        vy_rel = vy
        v_rel_mag = np.sqrt(vx_rel**2 + vy_rel**2)
        
        if v_rel_mag < 1e-6:
            Fx_drag = 0.0
            Fy_drag = 0.0
        else:
            Fx_drag = -0.5 * self.rho * Cd * A * v_rel_mag * vx_rel
            Fy_drag = -0.5 * self.rho * Cd * A * v_rel_mag * vy_rel
        
        # FIXED: Correct drag calculation
        ax = Fx_drag / m
        ay = -self.g + Fy_drag / m
        
        return [vx, ax, vy, ay]
    
    @staticmethod
    def _hit_ground(t, state, *args):
        """Event function: triggers when y = 0 (ground hit)."""
        return state[2]
    
    setattr(_hit_ground, 'terminal', True)
    setattr(_hit_ground, 'direction', -1)
    
    def simulate(
        self,
        wind_speed: float = 0.0,
        max_time: float = 1200.0,
    ) -> Tuple[Tuple[float, float], Tuple[float, float], Dict]:
        """
        Simulate complete two-stage rocket launch.
        
        Args:
            wind_speed: Horizontal wind speed (m/s)
            max_time: Maximum simulation time (s)
        
        Returns:
            payload_impact: (x, y) where payload hits ground
            stage1_impact: (x, y) where Stage 1 debris hits ground
            trajectory: Dictionary with full trajectory data
        
        Raises:
            RuntimeError: If integration fails
        """
        try:
            # Phase 1: Stage 1 burn
            state0 = np.array([0.0, 0.0, 0.0, 0.0])
            t_span1 = (0, self.stage1.burn_time)
            t_eval1 = np.linspace(*t_span1, 500)
            
            sol1 = solve_ivp(
                self._rocket_eom,
                t_span1,
                state0,
                args=(self.stage1, 0.0, wind_speed),
                t_eval=t_eval1,
                max_step=0.5,
                method="RK45",
            )
            
            if not sol1.success:
                raise RuntimeError(f"Stage 1 integration failed: {sol1.message}")
            
            # State at separation
            sep_state = sol1.y[:, -1]
            t_sep = self.stage1.burn_time
            
            # Stage 1 debris trajectory
            sol_debris = solve_ivp(
                self._ballistic_eom,
                (t_sep, t_sep + max_time),
                sep_state,
                args=(self.stage1.Cd, self.stage1.A, self.stage1.m_dry, wind_speed),
                events=self._hit_ground,
                max_step=0.5,
                method="RK45",
            )
            
            if not sol_debris.success and not sol_debris.t_events[0].size > 0:
                raise RuntimeError(f"Debris integration failed: {sol_debris.message}")
            
            if sol_debris.t_events[0].size > 0:
                stage1_impact = (sol_debris.y_events[0][0][0], 0.0)
            else:
                stage1_impact = (sol_debris.y[0, -1], 0.0)
            
            # Phase 2: Stage 2 burn
            t_span2 = (t_sep, t_sep + self.stage2.burn_time)
            t_eval2 = np.linspace(*t_span2, 500)
            
            sol2 = solve_ivp(
                self._rocket_eom,
                t_span2,
                sep_state,
                args=(self.stage2, t_sep, wind_speed),
                t_eval=t_eval2,
                max_step=0.5,
                method="RK45",
            )
            
            if not sol2.success:
                raise RuntimeError(f"Stage 2 integration failed: {sol2.message}")
            
            burnout_state = sol2.y[:, -1]
            t_burnout = t_sep + self.stage2.burn_time
            
            # Phase 3: Ballistic flight
            sol3 = solve_ivp(
                self._ballistic_eom,
                (t_burnout, t_burnout + max_time),
                burnout_state,
                args=(self.stage2.Cd, self.stage2.A, self.stage2.m_dry, wind_speed),
                events=self._hit_ground,
                max_step=0.5,
                method="RK45",
            )
            
            if not sol3.success and not sol3.t_events[0].size > 0:
                raise RuntimeError(f"Ballistic integration failed: {sol3.message}")
            
            if sol3.t_events[0].size > 0:
                payload_impact = (sol3.y_events[0][0][0], 0.0)
            else:
                payload_impact = (sol3.y[0, -1], 0.0)
            
            # Combine trajectories
            trajectory = {
                "stage1_burn": {"t": sol1.t, "x": sol1.y[0], "y": sol1.y[2]},
                "stage2_burn": {"t": sol2.t, "x": sol2.y[0], "y": sol2.y[2]},
                "ballistic": {"t": sol3.t, "x": sol3.y[0], "y": sol3.y[2]},
                "debris": {"t": sol_debris.t, "x": sol_debris.y[0], "y": sol_debris.y[2]},
            }
            
            return payload_impact, stage1_impact, trajectory
        
        except Exception as e:
            raise RuntimeError(f"Simulation failed: {str(e)}")
    
    def monte_carlo(
        self,
        n_samples: int = 10000,
        seed: int = 42,
        uncertainties: Optional[Dict] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if uncertainties is None:
            uncertainties = {
                "thrust_stage1_percent": 3.0,
                "thrust_stage2_percent": 3.0,
                "mass_stage1_percent": 1.0,
                "wind_speed_std": 5.0,
                "launch_angle_std_deg": 0.5,
            }
        
        rng = default_rng(seed)
        
        payload_impacts = np.zeros(n_samples)
        stage1_impacts = np.zeros(n_samples)
        
        for i in range(n_samples):
            try:
                # Sample uncertain parameters
                s1 = StageParameters(
                    m_total=rng.normal(
                        self.stage1.m_total,
                        self.stage1.m_total * uncertainties["mass_stage1_percent"] / 100.0
                    ),
                    m_propellant=self.stage1.m_propellant,
                    thrust=rng.normal(
                        self.stage1.thrust,
                        self.stage1.thrust * uncertainties["thrust_stage1_percent"] / 100.0
                    ),
                    burn_time=self.stage1.burn_time,
                    Cd=self.stage1.Cd,
                    A=self.stage1.A,
                )
                
                s2 = StageParameters(
                    m_total=self.stage2.m_total,
                    m_propellant=self.stage2.m_propellant,
                    thrust=rng.normal(
                        self.stage2.thrust,
                        self.stage2.thrust * uncertainties["thrust_stage2_percent"] / 100.0
                    ),
                    burn_time=self.stage2.burn_time,
                    Cd=self.stage2.Cd,
                    A=self.stage2.A,
                )
                
                wind = rng.normal(0.0, uncertainties["wind_speed_std"])
                theta_offset = rng.normal(
                    0.0,
                    np.radians(uncertainties["launch_angle_std_deg"])
                )
                theta_sample = self.theta + theta_offset
                
                # Clamp angle to valid range
                theta_sample = np.clip(theta_sample, 0, np.pi/2)
                
                # Create temporary simulator with sampled parameters
                sim = RocketSimulator(s1, s2, theta_sample, self.g, self.rho)
                payload, debris, _ = sim.simulate(wind_speed=wind)
                
                payload_impacts[i] = payload[0]
                stage1_impacts[i] = debris[0]
            
            except Exception:
                # Mark failed simulation as NaN
                payload_impacts[i] = np.nan
                stage1_impacts[i] = np.nan
            
            # Progress indicator
            if (i + 1) % 1000 == 0:
                print(f"  {i + 1}/{n_samples} simulations complete")
        
        return payload_impacts, stage1_impacts