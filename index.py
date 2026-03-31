import json
import numpy as np
import matplotlib.pyplot as plt
from rocket_simulator import RocketSimulator, StageParameters


# ============================================================
# CONFIGURATION LOADING
# ============================================================

def load_config(config_file: str = "config.json") -> dict:
    """Load simulation configuration from JSON file."""
    try:
        with open(config_file, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file '{config_file}' not found")
    except json.JSONDecodeError:
        raise ValueError(f"Invalid JSON in '{config_file}'")


# Visualization 

def plot_single_trajectory(simulator: RocketSimulator, wind_speed: float = 0.0):
    """Plot a single nominal trajectory."""
    print("Simulating nominal trajectory...")
    payload, debris, traj = simulator.simulate(wind_speed=wind_speed)
    
    fig, ax = plt.subplots(figsize=(12, 6), facecolor="black")
    ax.set_facecolor("black")
    
    # Combine all trajectory segments and filter y >= 0 (only above ground)
    x_all = np.concatenate([
        traj["stage1_burn"]["x"],
        traj["stage2_burn"]["x"],
        traj["ballistic"]["x"],
    ])
    y_all = np.concatenate([
        traj["stage1_burn"]["y"],
        traj["stage2_burn"]["y"],
        traj["ballistic"]["y"],
    ])
    
    # Filter: only plot above ground
    valid_idx = y_all >= 0
    x_plot = x_all[valid_idx]
    y_plot = y_all[valid_idx]
    
    # Debris also filtered
    debris_valid = traj["debris"]["y"] >= 0
    debris_x = traj["debris"]["x"][debris_valid]
    debris_y = traj["debris"]["y"][debris_valid]
    
    ax.plot(x_plot / 1000, y_plot / 1000,
            color="cyan", linewidth=1.5, label="Payload trajectory")
    ax.plot(debris_x / 1000, debris_y / 1000,
            color="red", linewidth=1, alpha=0.7, linestyle="--",
            label="Stage 1 debris")
    
    # Mark impact points
    ax.plot(payload[0] / 1000, 0, "o", color="cyan", markersize=8, label=f"Payload impact: {payload[0]/1000:.2f} km")
    ax.plot(debris[0] / 1000, 0, "x", color="red", markersize=10, markeredgewidth=2, label=f"Debris impact: {debris[0]/1000:.2f} km")
    
    ax.set_xlabel("Downrange (km)", color="white", fontsize=12)
    ax.set_ylabel("Altitude (km)", color="white", fontsize=12)
    ax.set_title("Nominal trajectory — two-stage rocket", color="white", fontsize=14)
    ax.tick_params(colors="white")
    ax.legend(facecolor="black", edgecolor="white", labelcolor="white", fontsize=9)
    ax.grid(alpha=0.15, color="white")
    
    plt.tight_layout()
    plt.savefig("nominal_trajectory.png", dpi=150, facecolor="black")
    print("Saved: nominal_trajectory.png")
    print("\nIMPACT POSITIONS:")
    print(f"  Payload: {payload[0]/1000:.2f} km downrange")
    print(f"  Debris:  {debris[0]/1000:.2f} km downrange")
    plt.show()


def plot_monte_carlo(payload_impacts: np.ndarray, stage1_impacts: np.ndarray):
    """Plot Monte Carlo results: scatter, histogram, statistics."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), facecolor="black")
    
    # --- Plot 1: Impact scatter ---
    ax1 = axes[0]
    ax1.set_facecolor("black")
    valid_p = payload_impacts[~np.isnan(payload_impacts)]
    valid_s = stage1_impacts[~np.isnan(stage1_impacts)]
    
    ax1.scatter(valid_p / 1000, np.zeros_like(valid_p),
                c="cyan", s=1, alpha=0.3, label="Payload")
    ax1.scatter(valid_s / 1000, np.ones_like(valid_s) * 0.1,
                c="red", s=1, alpha=0.3, label="Stage 1 debris")
    ax1.set_xlabel("Downrange (km)", color="white")
    ax1.set_title("Impact points", color="white", fontsize=13)
    ax1.tick_params(colors="white")
    ax1.legend(facecolor="black", edgecolor="white", labelcolor="white", fontsize=9)
    
    # --- Plot 2: Payload histogram ---
    ax2 = axes[1]
    ax2.set_facecolor("black")
    ax2.hist(valid_p / 1000, bins=80, color="cyan", alpha=0.7, edgecolor="none")
    
    mean_p = np.mean(valid_p) / 1000
    std_p = np.std(valid_p) / 1000
    ci95 = np.percentile(valid_p / 1000, [2.5, 97.5])
    
    ax2.axvline(mean_p, color="white", linestyle="-", linewidth=1, label=f"Mean: {mean_p:.1f} km")
    ax2.axvline(ci95[0], color="yellow", linestyle="--", linewidth=1, label=f"95% CI: {ci95[0]:.1f}–{ci95[1]:.1f} km")
    ax2.axvline(ci95[1], color="yellow", linestyle="--", linewidth=1)
    
    ax2.set_xlabel("Downrange (km)", color="white")
    ax2.set_ylabel("Count", color="white")
    ax2.set_title("Payload impact distribution", color="white", fontsize=13)
    ax2.tick_params(colors="white")
    ax2.legend(facecolor="black", edgecolor="white", labelcolor="white", fontsize=9)
    
    # --- Plot 3: Stage 1 debris histogram ---
    ax3 = axes[2]
    ax3.set_facecolor("black")
    ax3.hist(valid_s / 1000, bins=80, color="red", alpha=0.7, edgecolor="none")
    
    mean_s = np.mean(valid_s) / 1000
    std_s = np.std(valid_s) / 1000
    ci95_s = np.percentile(valid_s / 1000, [2.5, 97.5])
    
    ax3.axvline(mean_s, color="white", linestyle="-", linewidth=1, label=f"Mean: {mean_s:.1f} km")
    ax3.axvline(ci95_s[0], color="yellow", linestyle="--", linewidth=1, label=f"95% CI: {ci95_s[0]:.1f}–{ci95_s[1]:.1f} km")
    ax3.axvline(ci95_s[1], color="yellow", linestyle="--", linewidth=1)
    
    ax3.set_xlabel("Downrange (km)", color="white")
    ax3.set_ylabel("Count", color="white")
    ax3.set_title("Stage 1 debris distribution", color="white", fontsize=13)
    ax3.tick_params(colors="white")
    ax3.legend(facecolor="black", edgecolor="white", labelcolor="white", fontsize=9)
    
    plt.tight_layout()
    plt.savefig("monte_carlo_results.png", dpi=150, facecolor="black")
    print("Saved: monte_carlo_results.png")
    plt.show()
    
    # --- Print statistics ---
    print("\n" + "=" * 60)
    print("MONTE CARLO FLIGHT SAFETY ANALYSIS".center(60))
    print("=" * 60)
    print(f"\nTotal samples:    {len(valid_p)} (succeeded) / {len(payload_impacts)} (total)")
    print(f"Failed runs:      {np.sum(np.isnan(payload_impacts))}")
    
    print(f"\n{'PAYLOAD IMPACT DISTRIBUTION':^60}")
    print("-" * 60)
    print(f"  Mean:           {mean_p:>8.2f} km")
    print(f"  Std Dev:        {std_p:>8.2f} km")
    print(f"  95% CI:         {ci95[0]:>8.2f} – {ci95[1]:<8.2f} km")
    print(f"  Min:            {np.min(valid_p)/1000:>8.2f} km")
    print(f"  Max:            {np.max(valid_p)/1000:>8.2f} km")
    
    print(f"\n{'STAGE 1 DEBRIS DISTRIBUTION':^60}")
    print("-" * 60)
    print(f"  Mean:           {mean_s:>8.2f} km")
    print(f"  Std Dev:        {std_s:>8.2f} km")
    print(f"  95% CI:         {ci95_s[0]:>8.2f} – {ci95_s[1]:<8.2f} km")
    print(f"  Min:            {np.min(valid_s)/1000:>8.2f} km")
    print(f"  Max:            {np.max(valid_s)/1000:>8.2f} km")
    
    print(f"\n{'HAZARD AREAS (95% CONFIDENCE)':^60}")
    print("-" * 60)
    print(f"  Stage 1 debris: {ci95_s[0]:>8.2f} – {ci95_s[1]:<8.2f} km")
    print(f"  Payload:        {ci95[0]:>8.2f} – {ci95[1]:<8.2f} km")
    print("\n" + "=" * 60)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("ROCKET FLIGHT SIMULATION".center(60))
    print("=" * 60)
    
    # Load configuration
    print("\nLoading configuration...")
    config = load_config("config.json")
    
    # Build stage parameters
    s1_cfg = config["stage1"]
    s2_cfg = config["stage2"]
    sim_cfg = config["simulation"]
    env_cfg = config["environment"]
    mc_cfg = config["monte_carlo"]
    
    stage1 = StageParameters(
        m_total=s1_cfg["m_total"],
        m_propellant=s1_cfg["m_propellant"],
        thrust=s1_cfg["thrust"],
        burn_time=s1_cfg["burn_time"],
        Cd=s1_cfg["Cd"],
        A=s1_cfg["A"],
    )
    
    stage2 = StageParameters(
        m_total=s2_cfg["m_total"],
        m_propellant=s2_cfg["m_propellant"],
        thrust=s2_cfg["thrust"],
        burn_time=s2_cfg["burn_time"],
        Cd=s2_cfg["Cd"],
        A=s2_cfg["A"],
    )
    
    theta_rad = np.radians(sim_cfg["launch_angle_deg"])
    
    # Create simulator
    print("Initializing simulator...")
    simulator = RocketSimulator(
        stage1=stage1,
        stage2=stage2,
        launch_angle_rad=theta_rad,
        g=env_cfg["g"],
        rho=env_cfg["rho"],
    )
    
    # Plot nominal trajectory
    plot_single_trajectory(simulator, wind_speed=sim_cfg["wind_speed_nominal"])
    
    # Run Monte Carlo
    print(f"\nRunning Monte Carlo ({mc_cfg['n_samples']:,} samples)...")
    payload_impacts, stage1_impacts = simulator.monte_carlo(
        n_samples=mc_cfg["n_samples"],
        seed=mc_cfg["seed"],
        uncertainties=mc_cfg["uncertainties"],
    )
    
    # Plot results
    plot_monte_carlo(payload_impacts, stage1_impacts)

