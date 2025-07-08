import numpy as np
import cross_sections_data as csd
import monte_carlo_initialisations as mci
import calculating_geometric_properties as cgp
import energy_calculation as ec


def pair_production_simulation(position, energy, pztop, pzbottom, radius, compton_scattering_cross_sections, photoelectric_absorption_cross_sections):

    direction1 = mci.isotropic_direction_in_angle(np.pi)
    direction2 = -direction1
    energy1 = 0.511
    energy2 = 0.511
    position1 = position.copy()
    position2 = position.copy()

    compton_cross_section = csd.get_cross_section(energy1, energy, compton_scattering_cross_sections)
    photoelectric_cross_section = csd.get_cross_section(energy1, energy, photoelectric_absorption_cross_sections)
    total_cross_section = compton_cross_section + photoelectric_cross_section  

    summing_energy = 0

    while True:

        if energy1 < 0.001:
                break

        random_numbers = np.random.rand(2)

        lambd = -np.log(random_numbers[0]) / total_cross_section

        if cgp.goes_outside(position1, direction1, pztop, pzbottom, radius, lambd):
            break

        position1 += direction1 * lambd

        compton_scatter = random_numbers[1] < (compton_cross_section / total_cross_section)
        if compton_scatter:
            direction1, energy1, detector_energy = mci.compton_scatter_photon(energy1, direction1)

            summing_energy += detector_energy

        photoelectric_absorb = random_numbers[1] < ((photoelectric_cross_section + compton_cross_section) / total_cross_section) and not compton_scatter
        if photoelectric_absorb:

            summing_energy += energy1
            break


        compton_cross_section = csd.get_cross_section(energy1, energy, compton_scattering_cross_sections)
        photoelectric_cross_section = csd.get_cross_section(energy1, energy, photoelectric_absorption_cross_sections)
        total_cross_section = compton_cross_section + photoelectric_cross_section


    compton_cross_section = csd.get_cross_section(energy2, energy, compton_scattering_cross_sections)
    photoelectric_cross_section = csd.get_cross_section(energy2, energy, photoelectric_absorption_cross_sections)
    total_cross_section = compton_cross_section + photoelectric_cross_section  

    while True:

        if energy2 < 0.001:
            break
        
        random_numbers = np.random.rand(2)

        lambd = -np.log(random_numbers[0]) / total_cross_section

        if cgp.goes_outside(position2, direction2, pztop, pzbottom, radius, lambd):
            break

        position2 += direction2 * lambd

        compton_scatter = random_numbers[1] < (compton_cross_section / total_cross_section)
        if compton_scatter:
            direction2, energy2, detector_energy = mci.compton_scatter_photon(energy2, direction2)

            summing_energy += detector_energy

        photoelectric_absorb = random_numbers[1] < ((photoelectric_cross_section + compton_cross_section) / total_cross_section) and not compton_scatter
        if photoelectric_absorb:

            summing_energy += energy2
            break


        compton_cross_section = csd.get_cross_section(energy2, energy, compton_scattering_cross_sections)
        photoelectric_cross_section = csd.get_cross_section(energy2, energy, photoelectric_absorption_cross_sections)
        total_cross_section = compton_cross_section + photoelectric_cross_section

    return summing_energy


def simulate_transport(source_position, source_energy, energy, compton_scattering_cross_sections,
                        photoelectric_absorption_cross_sections, pair_production_cross_sections,
                        pztop, pzbottom, radius, FWHM, num_photons, alpha, detector_height, detector_radius):

    energy_accumulator = ec.EnergyHistogram(num_bins=1024, energy_min=0., energy_max = 1.1 * source_energy)

    reached_detector_num = 0


    for _ in range(num_photons):

        current_energy = source_energy

        direction = mci.isotropic_direction_in_cone(alpha, source_position, detector_height, detector_radius)

        if np.sqrt(source_position[0]**2 + source_position[1]**2) < radius and source_position[2] < pztop and source_position[2] > pzbottom:  
            position = source_position.copy() 
        else:
            position = cgp.intersect_cylinder_starting_points(source_position, direction, pztop, pzbottom, radius)

        summing_energy = 0

        if position is not None: #reached the detector
            reached_detector_num += 1


        while position is not None and current_energy > 0.001:

            compton = csd.get_cross_section(current_energy, energy, compton_scattering_cross_sections)
            photo = csd.get_cross_section(current_energy, energy, photoelectric_absorption_cross_sections)
            pair = csd.get_cross_section(current_energy, energy, pair_production_cross_sections) if current_energy > 1.022 else 0
            total = compton + photo + pair

            if total <= 0:
                break  

            lambd = -np.log(np.random.rand()) / total

            if cgp.goes_outside(position, direction, pztop, pzbottom, radius, lambd):
                if summing_energy > 0:
                    energy_accumulator.add([summing_energy * np.random.normal(1, FWHM/2.3548)])
                break

            position += direction * lambd


            rand = np.random.rand() * total
            if rand < compton:
                direction, current_energy, deposited = mci.compton_scatter_photon(current_energy, direction)
                summing_energy += deposited
                
            elif rand < compton + photo:
                summing_energy += current_energy
                break
                
            elif current_energy > 1.022:
                pair_energy = pair_production_simulation(position, energy, pztop, pzbottom, radius,
                                                    compton_scattering_cross_sections, 
                                                    photoelectric_absorption_cross_sections)
                summing_energy += current_energy - 1.022 + pair_energy
                break

        if summing_energy > 0:
            energy_accumulator.add([summing_energy * np.random.normal(1, FWHM / 2.3548)])

    E_int = reached_detector_num * source_energy

    return energy_accumulator, E_int


def record_gamma_spectrum(source_position, source_energy, detector_height, detector_radius, NaI_density, FWHM, num_particles, cross_sections_file_path=None):

    if cross_sections_file_path is None:
        raise ValueError("Cross sections file path must be provided.")

    energy, compton_scattering_cross_sections, photoelectric_absorption_cross_sections, pair_production_cross_sections = csd.read_cross_sections(cross_sections_file_path)

    energy = csd.make_energy_unique(energy)



    compton_scattering_cross_sections *= NaI_density
    photoelectric_absorption_cross_sections *= NaI_density
    pair_production_cross_sections *= NaI_density

    alpha = cgp.calc_angle_of_cone(source_position, detector_height, detector_radius)

    pztop = detector_height / 2
    pzbottom = -detector_height / 2

    energy_accumulator, E_int = simulate_transport(source_position, source_energy, energy, compton_scattering_cross_sections,
                        photoelectric_absorption_cross_sections, pair_production_cross_sections,
                        pztop, pzbottom, detector_radius, FWHM, num_particles, alpha, detector_height, detector_radius)
    
    E_det = ec.total_energy_in_histogram(energy_accumulator)

    if np.linalg.norm(source_position) < np.sqrt(detector_radius**2 + detector_height**2 / 4):
        E_tot = num_particles * source_energy
    else:
        E_tot = (2 / (1 - np.cos(alpha))) * num_particles * source_energy

    efficiency_tot = E_det / E_tot
    efficiency_int = E_det / E_int

    return energy_accumulator, efficiency_tot, efficiency_int