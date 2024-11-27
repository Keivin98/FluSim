# FluSim -- SEIRV Monte Carlo Simulation Library

This library provides a framework to simulate the spread of infectious diseases using an enhanced SEIRV (Susceptible, Exposed, Infectious, Recovered, Vaccinated, Dead) model. The simulations are run with Monte Carlo methods to account for stochastic variations, allowing users to evaluate the effects of interventions like masking, social distancing, and vaccination.

---

## Features

- **Modular Design**: Customize parameters for scenarios or use predefined ones.
- **Monte Carlo Simulations**: Perform multiple runs to compute confidence intervals for results.
- **Dynamic Infection Modeling**: Simulate realistic disease spread with proximity-based infection probabilities.
- **Vaccination Effects**: Model single-dose and two-dose vaccine efficacy with delays.
- **Flexible Output**: Generate detailed reports and visualize SEIRV dynamics with customizable plots.

---

# Simulation Results

Below is a grid showing the SEIRV simulation results for the predefined scenarios.

| **Baseline**                                           | **Masking and Social Distancing**                                                                |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------ |
| ![No Masking](plots/seirv_simulation_No%20Masking.png) | ![Masking and Social Distancing](plots/seirv_simulation_Masking%20and%20Social%20Distancing.png) |

| **1 Vaccine Dose**                                               | **2 Vaccine Doses**                                                |
| ---------------------------------------------------------------- | ------------------------------------------------------------------ |
| ![1 Vaccine Dose](plots/seirv_simulation_1%20Vaccine%20Dose.png) | ![2 Vaccine Doses](plots/seirv_simulation_2%20Vaccine%20Doses.png) |

---

## Installation

Clone the repository or copy the `main.py` file into your project.

```bash
git clone https://github.com/Keivin98/FluSim.git
cd FluSim
```

Install the required Python libraries:

```
pip install numpy matplotlib tqdm scipy
```

Usage

Running Predefined Scenarios

The library includes four predefined scenarios: 1. No Masking: No interventions. 2. Masking and Social Distancing: Masking reduces infection probability. 3. 1 Vaccine Dose: Vaccination with a single dose. 4. 2 Vaccine Doses: Masking combined with a two-dose vaccination strategy.

Run the simulation with:

```
python main.py
```

Results include:
• SEIRV plots saved as PNG files.
• Final reports for each scenario.

Customizing Simulations

Example: Custom Parameters

You can override default parameters by defining your own:

```
custom_params = {
    "population_size": 2000,
    "initial_infected": 50,
    "vaccination_rate": 0.02,
    "vaccination_delay": 14,
    "delay_between_doses": 21,
    "vaccine_efficacy_per_dose": [0.6, 0.3],
    "recovery_mean": 10,
    "recovery_sd": 2,
    "latent_period_mean": 2,
    "latent_period_sd": 1,
    "mortality_rate": 0.03,
    "contact_rate_lambda": 20,
    "masking_effectiveness": 0.4,
    "vaccination_start_step": 10,
    "total_steps": 200,
    "infection_radius": 0.075,
    "base_infection_prob": 0.3,
}
```

Run the simulation:

```
from seirv_simulation import run_monte_carlo_simulation, plot_seirv_curves, generate_report

runs = 50
per_state_counts_mean, per_state_counts_std, final_counts_mean, final_counts_std, report_data_aggregate = run_monte_carlo_simulation(
    runs=runs,
    params=custom_params,
    vaccination_rate=0.02,
    doses=2,
    delay_between_doses=21,
)

plot_seirv_curves(per_state_counts_mean, per_state_counts_std, "Custom Scenario")
generate_report(final_counts_mean, final_counts_std, report_data_aggregate, "Custom Scenario")
```

Reports

Each scenario generates a report summarizing:

Final Counts:
• Susceptible (S)
• Exposed (E)
• Infectious (I)
• Recovered (R)
• Deceased (D)
• Vaccinated (V)

Aggregate Data:
• Total Infected
• Total Deceased
• Total Recovered
• Total Vaccinated

Plots

SEIRV dynamics are visualized with:
• Customizable y-axis ticks: Based on powers of 2 (e.g., 5, 10, 20, 40, …).
• Confidence Intervals: Shaded regions show ±1 standard deviation.
• Logarithmic Scaling: Makes trends across different orders of magnitude more readable.

Plots are saved automatically as PNG files.

Key Parameters

| Parameter                   | Description                                                            | Default Value |
| --------------------------- | ---------------------------------------------------------------------- | ------------- |
| `population_size`           | Total number of people in the simulation.                              | 2000          |
| `initial_infected`          | Number of initially infected individuals.                              | 50            |
| `vaccination_rate`          | Daily proportion of the population vaccinated.                         | 0.02          |
| `vaccination_delay`         | Days before vaccine-induced immunity becomes effective.                | 14            |
| `delay_between_doses`       | Minimum time between vaccine doses.                                    | 21            |
| `vaccine_efficacy_per_dose` | Efficacy of each vaccine dose.                                         | [0.6, 0.3]    |
| `recovery_mean`             | Mean recovery time in days.                                            | 10            |
| `recovery_sd`               | Standard deviation of recovery time.                                   | 2             |
| `latent_period_mean`        | Mean latent period (days before becoming infectious).                  | 2             |
| `latent_period_sd`          | Standard deviation of the latent period.                               | 1             |
| `mortality_rate`            | Probability of death for an infected individual.                       | 0.03          |
| `contact_rate_lambda`       | Average number of daily contacts per person.                           | 20            |
| `masking_effectiveness`     | Proportionate reduction in transmission probability due to masking.    | 0.4           |
| `vaccination_start_step`    | Simulation step when vaccination begins.                               | 10            |
| `total_steps`               | Total number of days in the simulation.                                | 200           |
| `infection_radius`          | Distance within which an infectious individual can spread the disease. | 0.075         |
| `base_infection_prob`       | Probability of transmission per contact (before adjustments).          | 0.3           |
