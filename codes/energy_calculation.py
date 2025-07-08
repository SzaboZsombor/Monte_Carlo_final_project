import numpy as np
import matplotlib.pyplot as plt


class EnergyHistogram:
    def __init__(self, num_bins=1024, energy_min=0.5, energy_max=4.2):
        self.num_bins = num_bins
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.bin_edges = np.linspace(energy_min, energy_max, num_bins + 1)
        self.hist = np.zeros(num_bins, dtype=int)
        self.energies = []

    def add(self, energies):
        counts, _ = np.histogram(energies, bins=self.bin_edges)
        self.hist += counts
        self.energies.extend(energies)

    def get_histogram(self):
        return self.hist, self.bin_edges
    
    def get_all_energies(self):
        return np.array(self.energies)
    

def total_energy_in_histogram(energy_histogram):
    
    energies = energy_histogram.get_all_energies()

    return np.sum(energies)


def plot_energy_accumulator(energy_accumulator, title='Gamma Spectrum Monte Carlo Simulation', E_gamma = 0., save_path=None):

    m_e = 0.511

    compton_edge = E_gamma * (1 - 1 / (1 + 2 * E_gamma / m_e))

    single_escape = E_gamma - 0.511 if E_gamma > 1.022 else None

    double_escape = E_gamma - 1.022 if E_gamma > 1.022 else None

    hist, bin_edges = energy_accumulator.get_histogram()

    plt.figure(figsize=(10, 6))
    plt.hist(
        bin_edges[:-1],  
        bins=bin_edges,          
        weights=hist,            
        edgecolor='black',      
        linewidth=0.5,           
        alpha=0.7,               
        log=True                 
    )
    plt.title(title)
    plt.ylim(0.5, 1e4)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Counts [-]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.axvline(compton_edge, color='orange', linestyle='--', linewidth=2, label='Compton edge')
    if single_escape is not None:
        plt.axvline(single_escape, color='green', linestyle='--', linewidth=2, label='Pair production, single escape')
    if double_escape is not None:
        plt.axvline(double_escape, color='red', linestyle='--', linewidth=2, label='Pair production, double escape')
    plt.legend()
    if save_path:
        plt.savefig(save_path)
    plt.show()