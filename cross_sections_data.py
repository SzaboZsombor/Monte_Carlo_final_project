import numpy as np
import scipy.interpolate


def read_cross_sections(filepath):
    energy = []
    incoherent_scattering_cross_section = []
    photoelectric_absorption_cross_section = []
    nuclear_production_cross_section = []
    electron_production_cross_section = []
    with open(filepath, 'r') as file:
        i = 0
        for line in file:
            if i < 3:
                i += 1
            else:
                parts = line.split()
                if len(parts) == 5:
                    energy.append(float(parts[0]))
                    incoherent_scattering_cross_section.append(float(parts[1]))
                    photoelectric_absorption_cross_section.append(float(parts[2]))
                    nuclear_production_cross_section.append(float(parts[3]))
                    electron_production_cross_section.append(float(parts[4]))

    incoherent_scattering_cross_section = np.array(incoherent_scattering_cross_section)
    photoelectric_absorption_cross_section = np.array(photoelectric_absorption_cross_section)
    nuclear_production_cross_section = np.array(nuclear_production_cross_section)
    electron_production_cross_section = np.array(electron_production_cross_section)
    energy = np.array(energy)

    pair_production_cross_section = nuclear_production_cross_section + electron_production_cross_section

    return energy, incoherent_scattering_cross_section, photoelectric_absorption_cross_section, pair_production_cross_section


def make_energy_unique(energy, epsilon=1e-9):
    energy = energy.copy()
    for i in range(1, len(energy)):
        while energy[i] <= energy[i - 1]:
            energy[i] = energy[i - 1] + epsilon
    return energy


def get_cross_section(current_energy, energy, cross_section):
    energy = np.array(energy)
    cross_section = np.array(cross_section)
    
    if len(energy) != len(cross_section):
        raise ValueError("Energy and cross section arrays must have the same length.")
    
    interpolated_cross_section = scipy.interpolate.interp1d(energy, cross_section, bounds_error=False, fill_value=0.0)
    
    return interpolated_cross_section(current_energy)

