from scipy.stats import f_oneway
import numpy as np
import random
from scipy.stats import poisson, norm, bernoulli
from tqdm import tqdm
import matplotlib.pyplot as plt

# Constants (can be overridden in functions)
DEFAULT_CONSTANTS = {
    "population_size": 10000,
    "initial_infected": 50,
    "vaccination_rate": 0.0075,
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
    "vaccination_start_step": 30,
    "total_steps": 90,
    "infection_radius": 0.001,
    "base_infection_prob": 0.4,
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
        self.state = SUSCEPTIBLE
        self.vaccinated = False
        self.doses_received = 0
        self.vaccination_times = []
        self.vaccine_efficacy = 0.0
        self.infected_time = 0
        self.latent_period = max(1, int(norm.rvs(params["latent_period_mean"], params["latent_period_sd"])))
        self.recovery_time = self.latent_period + max(1, int(norm.rvs(params["recovery_mean"], params["recovery_sd"])))
        self.was_infected = False
        self.recovered_step = None

    def update_vaccination_status(self, current_step, params):
        """Update vaccine efficacy and vaccination status based on the step."""
        for idx, dose_time in enumerate(self.vaccination_times):
            # Check if the dose has taken effect
            if current_step - dose_time >= params["vaccination_delay"] and idx >= self.doses_received:
                self.doses_received += 1
                if idx < len(params["vaccine_efficacy_per_dose"]):
                    dose_efficacy = params["vaccine_efficacy_per_dose"][idx]
                else:
                    dose_efficacy = 0.0
                self.vaccine_efficacy += dose_efficacy
                self.vaccine_efficacy = min(self.vaccine_efficacy, 1.0)

                # Mark as vaccinated if the first dose is applied
                if not self.vaccinated:
                    self.vaccinated = True

def vaccinate_population(step, population, params, doses, delay_between_doses):
    """Select eligible people for vaccination and schedule doses."""
    num_to_vaccinate = int(params["vaccination_rate"] * params["population_size"])
    eligible_people = []
    for person in population:
        # Only vaccinate susceptible people who haven't exceeded the maximum doses
        if (
            person.state == SUSCEPTIBLE 
            and len(person.vaccination_times) < doses 
            and (not person.vaccination_times or (step - person.vaccination_times[-1]) >= delay_between_doses)
        ):
            eligible_people.append(person)

    if eligible_people:
        to_vaccinate = random.sample(eligible_people, min(num_to_vaccinate, len(eligible_people)))
        for person in to_vaccinate:
            person.vaccination_times.append(step)

def run_simulation(params=DEFAULT_CONSTANTS, vaccination_rate=None, masking_effectiveness=None,
                   infection_radius=None, doses=1, delay_between_doses=0, collect_data=True):
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
        person.state = EXPOSED
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
                        base_prob = params["base_infection_prob"]
                        vaccine_efficacy_susceptible = other_person.vaccine_efficacy
                        masking_effectiveness = params["masking_effectiveness"]
                        infection_prob = base_prob * (1 - vaccine_efficacy_susceptible) * (1 - masking_effectiveness)
                        random_prob = random.random()
                        if random_prob < infection_prob:
                            # print(random_prob, infection_prob)
                            other_person.state = EXPOSED
                            other_person.infected_time = 0
                            other_person.was_infected = True

        # Update health status
        for person in population:
            if person.state == EXPOSED:
                person.infected_time += 1
                if person.infected_time >= person.latent_period:
                    person.state = INFECTIOUS
                    person.infected_time = 0
            elif person.state == INFECTIOUS:
                person.infected_time += 1
                if person.infected_time >= (person.recovery_time - person.latent_period):
                    if bernoulli.rvs(params["mortality_rate"]):
                        person.state = DEAD
                    else:
                        person.state = RECOVERED
                        person.recovered_step = step

        # Collect per-state counts
        if collect_data:
            counts = {
                SUSCEPTIBLE: 0,
                EXPOSED: 0,
                INFECTIOUS: 0,
                RECOVERED: 0,
                DEAD: 0,
                'V': 0
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
    total_steps = params.get("total_steps", 90)
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
    total_steps = params.get("total_steps", 90)
    states = [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED, DEAD, 'V']
    time = np.arange(total_steps)
    
    plt.figure(figsize=(12, 8))
    
    for state in states:
        mean = per_state_counts_mean[state]
        std = per_state_counts_std[state]
        plt.plot(time, mean, label=state, linewidth=2)
        plt.fill_between(time, mean - std, mean + std, alpha=0.2)
    
    plt.xlabel('Time (days)', fontsize=20)
    plt.ylabel('Number of People', fontsize=20)
    plt.title(f'SEIRV Model Simulation - {scenario_name}', fontsize=32)
    
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.legend(fontsize=22)
    
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    
    if save:
        plt.savefig(f'plots/seirv_simulation_{scenario_name}.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()
    
    plt.close()

def generate_report(final_counts_mean, final_counts_std, report_data_aggregate, scenario_name):
    print(f"\n=== Report for {scenario_name} ===")

    print("\nAggregate Report Data:")
    for key, value in report_data_aggregate.items():
        print(f"{key}: {value['mean']:.2f} ± {value['std']:.2f}")
    print("="*30)

def perform_pairwise_anova(data, metric, pair_labels):
    print(f"\nPerforming Pairwise ANOVA for {metric}:")
    results = {}
    for (scenario1, scenario2) in pair_labels:
        group1 = data[scenario1][metric]
        group2 = data[scenario2][metric]
        
        # Perform ANOVA between the two groups
        f_stat, p_value = f_oneway(group1, group2)
        
        # Store and print the results
        results[f"{scenario1} vs. {scenario2}"] = (f_stat, p_value)
        print(f"{scenario1} vs. {scenario2}: F-statistic = {f_stat:.2f}, p-value = {p_value:.3e}")
    return results

if __name__ == "__main__":
    runs = 50  # Number of simulation runs for Monte Carlo
    
    scenarios = {
        "No Masking": {
            "masking_effectiveness": 0.0,
            "vaccination_rate": 0.0,
            "doses": 0,
            "delay_between_doses": 0,
        },
        "Masking and Social Distancing": {
            "masking_effectiveness": 0.45,
            "vaccination_rate": 0.0,
            "doses": 0,
            "delay_between_doses": 0,
        },
        "1 Vaccine Dose": {
            "masking_effectiveness": 0.45,
            "vaccination_rate": 0.02,
            "doses": 1,
            "delay_between_doses": 7,
        },
        "2 Vaccine Doses": {
            "masking_effectiveness": 0.45,
            "vaccination_rate": 0.02,
            "doses": 2,
            "delay_between_doses": 7,
        },
    }
    results = {}
    for scenario_name, scenario_params in scenarios.items():
        print(f"\nRunning Scenario: {scenario_name}")
        
        params = {
            "population_size": 10000,
            "total_steps": 200,
            "vaccination_start_step": 15,
            "vaccine_efficacy_per_dose": [0.4, 0.4],
            "vaccination_delay": 0,
            "delay_between_doses": scenario_params["delay_between_doses"],
            "contact_rate_lambda": 20,
            "recovery_mean": 10,
            "masking_effectiveness": scenario_params["masking_effectiveness"],
            "infection_radius": 0.1,
            "base_infection_prob": 0.4,
        }

        per_state_counts_mean, per_state_counts_std, final_counts_mean, final_counts_std, report_data_aggregate = run_monte_carlo_simulation(
            runs=runs,
            params=params,
            vaccination_rate=scenario_params["vaccination_rate"],
            doses=scenario_params["doses"],
            delay_between_doses=scenario_params["delay_between_doses"],
        )

        results[scenario_name] = {
            "infected": np.random.normal(final_counts_mean["I"], final_counts_std["I"], runs),
            "deceased": np.random.normal(final_counts_mean["D"], final_counts_std["D"], runs),
        }

        # Plot SEIRV curves
        plot_seirv_curves(per_state_counts_mean, per_state_counts_std, scenario_name, save=True, params=params)

        # Generate report
        generate_report(final_counts_mean, final_counts_std, report_data_aggregate, scenario_name)
    
    pair_labels = [
        ("No Masking", "Masking and Social Distancing"),
        ("Masking and Social Distancing", "1 Vaccine Dose"),
        ("1 Vaccine Dose", "2 Vaccine Doses"),
    ]

    # Perform ANOVA for infections and deaths
    infections_results = perform_pairwise_anova(results, "infected", pair_labels)
    deaths_results = perform_pairwise_anova(results, "deceased", pair_labels)