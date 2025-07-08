import numpy as np
import matplotlib.pyplot as plt
from energy_calculation import plot_energy_accumulator
from transport_simulation import record_gamma_spectrum

def plot_geometries(positions, detector_height, detector_radius, save_path=None):

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    z = np.linspace(-detector_height/2, detector_height/2, 50)
    theta = np.linspace(0, 2 * np.pi, 50)
    theta_grid, z_grid = np.meshgrid(theta, z)
    x_grid = detector_radius * np.cos(theta_grid)
    y_grid = detector_radius * np.sin(theta_grid)
    ax.plot_surface(x_grid, y_grid, z_grid, alpha=0.2, color='blue', linewidth=0)

    positions = np.array(positions)
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], color='red', s=50, label='Source positions')

    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_zlabel('Z [cm]')
    ax.set_title('Detector Geometry and Source Positions')
    ax.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()


def main():
    source_A = np.array([3.0, -3.0, 2.0])
    source_B = np.array([4.0, 4.0, 0.0])
    energy_A = 0.6617  # MeV
    energy_B = 1.3325  # MeV
    radius_A = 2.5
    radius_B = 3.0
    height_A = 3.0
    height_B = 5.0
    density = 3.67
    FWHM_A = 6.0 / 1000  # MeV
    FWHM_B = 8.0 / 1000  # MeV
    num_particles = 500000

    cross_sections_file_path = r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\cross_sections_data.txt'
    energy_acc_A, eff_tot_A, eff_int_A = record_gamma_spectrum(source_A, energy_A, height_A, radius_A, density, FWHM_A, num_particles, cross_sections_file_path=cross_sections_file_path)
    plot_energy_accumulator(energy_acc_A, f'Gamma Spectrum ({energy_A} MeV)', energy_A, save_path=r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\energy_spectrum_A.png')
    print(f"Spectrum A ({energy_A} MeV): η_tot =", eff_tot_A, "η_int =", eff_int_A)

    energy_acc_B, eff_tot_B, eff_int_B = record_gamma_spectrum(source_B, energy_B, height_B, radius_B, density, FWHM_B, num_particles, cross_sections_file_path=cross_sections_file_path)
    plot_energy_accumulator(energy_acc_B, f'Gamma Spectrum ({energy_B} MeV)', energy_B, save_path=r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\energy_spectrum_B.png')
    print(f"Spectrum B ({energy_B} MeV): η_tot =", eff_tot_B, "η_int =", eff_int_B)

    positions = np.linspace([1.0, 3.5, 2.0], [-4.0, -1.5, 2.0], 11)
    effs_tot = []
    effs_int = []

    for pos in positions:
        _, eff_tot, eff_int = record_gamma_spectrum(pos, energy_A, height_A, radius_A, density, FWHM_A, num_particles, cross_sections_file_path=cross_sections_file_path)
        effs_tot.append(eff_tot)
        effs_int.append(eff_int)


    plot_geometries(positions, height_A, radius_A, save_path=r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\detector_geometry.png')

    plt.figure()
    plt.plot(range(1, 12), effs_tot, label=r'$η_tot$')
    plt.plot(range(1, 12), effs_int, label=r'$η_int$')
    plt.xlabel('Position number')
    plt.ylabel('Efficiency')
    plt.title('Efficiencies vs. Source Position (Spectrum A)')
    plt.legend()
    plt.grid(True)
    plt.savefig(r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\efficiencies_vs_position.png')
    plt.show()


    energies = np.linspace(0.4, 4.0, 10)
    effs_tot_B = []
    effs_int_B = []

    for E in energies:
        _, eff_tot, eff_int = record_gamma_spectrum(source_B, E, height_B, radius_B, density, FWHM_B, num_particles, cross_sections_file_path=cross_sections_file_path)
        effs_tot_B.append(eff_tot)
        effs_int_B.append(eff_int)

    plt.figure()
    plt.plot(energies, effs_tot_B, label=r'$η_tot$')
    plt.plot(energies, effs_int_B, label=r'$η_int$')
    plt.xlabel('Source energy (MeV)')
    plt.ylabel('Efficiency')
    plt.title('Efficiencies vs. Source Energy (Spectrum B)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.savefig(r'D:\Egyetem\Monte Carlo\Monte_Carlo_final_project\plots_and_data\efficiencies_vs_energy.png')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()