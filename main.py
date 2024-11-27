import numpy as np
import random
from scipy.stats import poisson, norm, bernoulli, beta
from tqdm import tqdm
import matplotlib.pyplot as plt

# Constants (can be overridden in functions)
DEFAULT_CONSTANTS = {
    "population_size": 2000,  # Reduced population for quicker demonstration
    "initial_infected": 50,
    "vaccination_rate": 0.005,  # Default vaccination rate
    "vaccination_delay": 14,
    "delay_between_doses": 21,
    "vaccine_efficacy_per_dose": [0.6, 0.3],  # Efficacy after each dose
    "recovery_mean": 10,
    "recovery_sd": 2,
    "latent_period_mean": 2,
    "latent_period_sd": 1,
    "mortality_rate": 0.03,
    "contact_rate_lambda": 20,  # Reduced contact rate
    "masking_effectiveness": 0.4,  # Default masking effectiveness
    "vaccination_start_step": 50,
    "total_steps": 200,        # Adjusted total steps
    "infection_radius": 0.075,  # Adjusted infection radius
    "base_infection_prob": 0.85,  # Reduced base infection probability
}

# States
SUSCEPTIBLE = 'S'
EXPOSED = 'E'
INFECTIOUS = 'I'
RECOVERED = 'R'
DEAD = 'D'

class Person:
    def __init__(self, x, y, params):
        self.x = x
        self.y = y
        self.state = SUSCEPTIBLE  # 'S', 'E', 'I', 'R', 'D'
        self.vaccinated = False
        self.doses_received = 0
        self.vaccination_times = []  # Times when doses were received
        self.vaccine_efficacy = 0.0  # Total efficacy
        self.infected_time = 0
        self.latent_period = max(1, int(norm.rvs(params["latent_period_mean"], params["latent_period_sd"])))
        self.recovery_time = self.latent_period + max(1, int(norm.rvs(params["recovery_mean"], params["recovery_sd"])))
        self.was_infected = False
        self.recovered_step = None  # To track immunity duration if needed

    def update_vaccination_status(self, current_step, params):
        # Check for each dose if the vaccination delay has passed
        for idx, dose_time in enumerate(self.vaccination_times):
            if current_step - dose_time == params["vaccination_delay"] and idx >= self.doses_received:
                self.doses_received += 1
                if idx < len(params["vaccine_efficacy_per_dose"]):
                    dose_efficacy = params["vaccine_efficacy_per_dose"][idx]
                else:
                    dose_efficacy = 0.0  # No efficacy if doses exceed defined efficacies
                self.vaccine_efficacy += dose_efficacy
                self.vaccine_efficacy = min(self.vaccine_efficacy, 1.0)
                self.vaccinated = True
                # Optionally, handle second dose differently if needed

    def get_infection_probability(self, params):
        # Adjust infection probability based on vaccine efficacy
        adjusted_prob = params["base_infection_prob"] * (1 - self.vaccine_efficacy)
        return adjusted_prob

def vaccinate_population(step, population, params, doses, delay_between_doses):
    num_to_vaccinate = int(params["vaccination_rate"] * params["population_size"])
    eligible_people = []
    for person in population:
        if person.state == SUSCEPTIBLE and person.doses_received < doses:
            if not person.vaccination_times or (step - person.vaccination_times[-1]) >= delay_between_doses:
                eligible_people.append(person)
    if eligible_people:
        to_vaccinate = random.sample(eligible_people, min(num_to_vaccinate, len(eligible_people)))
        for person in to_vaccinate:
            person.vaccination_times.append(step)

def categorize_positions(population):
    positions = {
        SUSCEPTIBLE: [],
        EXPOSED: [],
        INFECTIOUS: [],
        RECOVERED: [],
        DEAD: [],
    }
    for person in population:
        positions[person.state].append((person.x, person.y))
    return positions

def run_simulation(params=DEFAULT_CONSTANTS, vaccination_rate=None, masking_effectiveness=None,
                   infection_radius=None, doses=1, delay_between_doses=0, collect_data=True):
    # Allow overriding specific constants
    params = {**DEFAULT_CONSTANTS, **params}
    if vaccination_rate is not None:
        params["vaccination_rate"] = vaccination_rate
    if masking_effectiveness is not None:
        params["masking_effectiveness"] = masking_effectiveness
    if infection_radius is not None:
        params["infection_radius"] = infection_radius

    # Initialize population
    population = [Person(np.random.random(), np.random.random(), params) for _ in range(params["population_size"])]
    for i in range(params["initial_infected"]):
        person = population[i]
        person.state = EXPOSED  # Start in Exposed state
        person.was_infected = True

    if collect_data:
        per_state_counts = []

    for step in range(params["total_steps"]):
        # Update vaccination status
        for person in population:
            person.update_vaccination_status(step, params)

        # Vaccinate population
        if step >= params["vaccination_start_step"] and params["vaccination_rate"] > 0:
            vaccinate_population(step, population, params, doses, delay_between_doses)

        # Spread infection
        infectious_people = [p for p in population if p.state == INFECTIOUS]
        for person in infectious_people:
            contacts = poisson.rvs(params["contact_rate_lambda"])
            for _ in range(contacts):
                other_person = random.choice(population)
                if other_person.state == SUSCEPTIBLE:
                    distance = np.sqrt((person.x - other_person.x) ** 2 + (person.y - other_person.y) ** 2)
                    if distance < params["infection_radius"]:
                        infection_prob = person.get_infection_probability(params)
                        infection_prob *= (1 - params["masking_effectiveness"])
                        if random.random() < infection_prob:
                            other_person.state = EXPOSED
                            other_person.infected_time = 0
                            other_person.was_infected = True

        # Update health status
        for person in population:
            if person.state == EXPOSED:
                person.infected_time += 1
                if person.infected_time >= person.latent_period:
                    person.state = INFECTIOUS
                    person.infected_time = 0  # Reset for infectious period
            elif person.state == INFECTIOUS:
                person.infected_time += 1
                if person.infected_time >= (person.recovery_time - person.latent_period):
                    if bernoulli.rvs(params["mortality_rate"]):
                        person.state = DEAD
                    else:
                        person.state = RECOVERED
                        person.recovered_step = step  # Track when recovered

        # Collect per-state counts
        if collect_data:
            counts = {
                SUSCEPTIBLE: 0,
                EXPOSED: 0,
                INFECTIOUS: 0,
                RECOVERED: 0,
                DEAD: 0,
                'V': 0  # Vaccinated count
            }
            for person in population:
                counts[person.state] += 1
                if person.vaccinated:
                    counts['V'] += 1
            per_state_counts.append(counts)

    # Prepare report data
    report_data = {
        'total_infected': sum(1 for p in population if p.was_infected),
        'total_deceased': sum(1 for p in population if p.state == DEAD),
        'total_recovered': sum(1 for p in population if p.state == RECOVERED),
        'total_vaccinated': sum(1 for p in population if p.vaccinated),
    }

    if collect_data:
        return per_state_counts, report_data
    else:
        return report_data

def run_monte_carlo_simulation(runs, params=DEFAULT_CONSTANTS, **kwargs):
    total_steps = params.get("total_steps", 365)
    states = [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED, DEAD, 'V']
    per_state_counts_runs = {state: np.zeros((runs, total_steps)) for state in states}
    final_counts_runs = {state: np.zeros(runs) for state in states}
    report_data_runs = []

    for i in tqdm(range(runs), desc="Monte Carlo Simulation"):
        per_state_counts, report_data = run_simulation(params=params, collect_data=True, **kwargs)
        report_data_runs.append(report_data)
        for t in range(total_steps):
            counts = per_state_counts[t]
            for state in states:
                per_state_counts_runs[state][i, t] = counts.get(state, 0)
        final_counts = per_state_counts[-1]
        for state in states:
            final_counts_runs[state][i] = final_counts.get(state, 0)

    # Compute averages and standard deviations
    per_state_counts_mean = {state: np.mean(per_state_counts_runs[state], axis=0) for state in states}
    per_state_counts_std = {state: np.std(per_state_counts_runs[state], axis=0) for state in states}

    # Compute final counts statistics
    final_counts_mean = {state: np.mean(final_counts_runs[state]) for state in states}
    final_counts_std = {state: np.std(final_counts_runs[state]) for state in states}

    # Aggregate report data
    report_data_aggregate = {}
    for key in report_data_runs[0].keys():
        values = [report[key] for report in report_data_runs]
        report_data_aggregate[key] = {
            'mean': np.mean(values),
            'std': np.std(values)
        }

    return per_state_counts_mean, per_state_counts_std, final_counts_mean, final_counts_std, report_data_aggregate

def plot_seirv_curves(per_state_counts_mean, per_state_counts_std, scenario_name, save=False, params=DEFAULT_CONSTANTS):
    total_steps = params.get("total_steps", 365)
    states = [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED, DEAD, 'V']
    time = np.arange(total_steps)
    plt.figure(figsize=(12, 8))
    for state in states:
        mean = per_state_counts_mean[state]
        std = per_state_counts_std[state]
        plt.plot(time, mean, label=state)
        plt.fill_between(time, mean - std, mean + std, alpha=0.2)
    plt.xlabel('Time (days)')
    plt.ylabel('Number of People')
    plt.title(f'SEIRV Model Simulation - {scenario_name}')
    plt.legend()
    if save:
        plt.savefig(f'seirv_simulation_{scenario_name}.png')
    else: 
        plt.show()
    plt.close()

def generate_report(final_counts_mean, final_counts_std, report_data_aggregate, scenario_name):
    print(f"\n=== Report for {scenario_name} ===")
    print("Final Counts (Mean ± Std):")
    for state in [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED, DEAD, 'V']:
        mean = final_counts_mean[state]
        std = final_counts_std[state]
        print(f"{state}: {mean:.2f} ± {std:.2f}")
    print("\nAggregate Report Data:")
    for key, value in report_data_aggregate.items():
        print(f"{key}: {value['mean']:.2f} ± {value['std']:.2f}")
    print("="*30)

# Define the four scenarios
scenarios = {
    "No Masking": {
        "masking_effectiveness": 0.0,
        "vaccination_rate": 0.0,
        "doses": 0,  # No vaccination
        "delay_between_doses": 0,
    },
    "Masking and Social Distancing": {
        "masking_effectiveness": 0.4,
        "vaccination_rate": 0.0,
        "doses": 0,  # No vaccination
        "delay_between_doses": 0,
    },
    "1 Vaccine Dose": {
        "masking_effectiveness": 0.0,
        "vaccination_rate": 0.005,
        "doses": 1,
        "delay_between_doses": 0,  # Not applicable for 1 dose
    },
    "2 Vaccine Doses": {
        "masking_effectiveness": 0.4,
        "vaccination_rate": 0.005,
        "doses": 2,
        "delay_between_doses": 21,
    },
}

# Example Usage: Run all scenarios
if __name__ == "__main__":
    runs = 50  # Number of simulation runs for Monte Carlo

    for scenario_name, scenario_params in scenarios.items():
        print(f"\nRunning Scenario: {scenario_name}")
        
        # Prepare parameters for simulation
        params = {
            "population_size": 2000,
            "total_steps": 200,
            "vaccination_start_step": 50,
            "vaccine_efficacy_per_dose": [0.6, 0.3],
            "vaccination_delay": 14,
            "delay_between_doses": scenario_params["delay_between_doses"],
            "contact_rate_lambda": 20,
            "recovery_mean": 10,
            "masking_effectiveness": scenario_params["masking_effectiveness"],
            "infection_radius": 0.075,
            "base_infection_prob": 0.85,
        }

        # Run Monte Carlo simulation for the current scenario
        per_state_counts_mean, per_state_counts_std, final_counts_mean, final_counts_std, report_data_aggregate = run_monte_carlo_simulation(
            runs=runs,
            params=params,
            vaccination_rate=scenario_params["vaccination_rate"],
            doses=scenario_params["doses"],
            delay_between_doses=scenario_params["delay_between_doses"],
        )

        # Plot SEIRV curves
        plot_seirv_curves(per_state_counts_mean, per_state_counts_std, scenario_name, save=True, params=params)

        # Generate report
        generate_report(final_counts_mean, final_counts_std, report_data_aggregate, scenario_name)