from parameters2 import *
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import warnings
import os
import pandas

style.use("bmh")
warnings.simplefilter("ignore")

def efficiency(T_cell, T_ref=298, eta_ref=0.156, beta_ref=0.0045):

    # Return efficiency
    return eta_ref * (1 - beta_ref * (T_cell - T_ref))


def plot(t, x1, x2):

    fig1, ax = plt.subplots(figsize=(10, 10))

    # ax1 = fig1.add_subplot(111)
    line1 = ax.plot(t, x1, 'o-', c="black", label="Temperature")
    plt.yticks(fontsize=18)
    plt.ylabel("Temperature (K)", fontsize=24)
    plt.xlabel("Time (hr)", fontsize=24)
    x_ticks = range(0, int(t[-1] + 1), 60 * 60)
    plt.xticks(ticks=x_ticks, labels=[str(int(x / 3600)) for x in x_ticks], fontsize=18)
    ax.legend(frameon=False, loc="upper right")


    ax2 = fig1.add_subplot(111, sharex=ax, frameon=False)
    line2 = ax2.plot(t, x2, 'xr-', c="black", label="Efficiency")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    plt.yticks(fontsize=18)
    plt.ylabel("Efficiency (%)", fontsize=24)
    plt.xticks(ticks=None, labels=None, fontsize=18)
    ax2.legend(frameon=False, loc="upper left")

    plt.show()


def simulation(plotting, PCM_thickness, PCM_density, PCM_latent_heat, PCM_conductivity, T_melt):

    PCM_sol_density = PCM_density * 1.05
    PCM_liq_density = PCM_density * 0.95
    PCM_sol_conductivity = PCM_conductivity * 1.15  # [W/(m*K)]
    PCM_liq_conductivity = PCM_conductivity * 0.85  # [W/(m*K)]
    PCM_sol_heat_cap = 2900  # [J/(kg*K)]
    PCM_liq_heat_cap = 2100  # [J/(kg*K)]
    shape_factor = 2
    lower_boundary = 0.001 * PCM_thickness  # [m]
    upper_boundary = 0.999 * PCM_thickness  # [m]

    # Define time domain and initial states
    t_real = np.linspace(t_start, t_end, int(1 + (t_end - t_start) / (dt)))
    t = np.linspace(t_start, t_end_ints, int(1 + (t_end_ints - t_start) / (dt_ints)))

    T_cell = np.zeros(len(t))
    cell_efficiency = np.zeros(len(t))
    total_power = 0

    T_cell[0] = T_air
    T_bottom = T_air
    T_front = T_air
    x_front = 0.01 * PCM_thickness
    cell_efficiency[0] = efficiency(T_cell[0])

    all_values = dict(zip(t, zip(T_cell, cell_efficiency)))


    for n in t:

        T_cell, cell_efficiency = all_values[n]

        #Safety feature
        if T_cell < 290 or T_cell > 340:
            print("Something went wrong. Skipping to next set of parameters.")
            return None, None, None

        # Calculate the energy absorbed by the cell and converted to heat
        total_power += cell_absorptance * cell_efficiency * solar_irridance
        cell_absorbed = cell_absorptance * ((1 - cell_efficiency) * solar_irridance)

        # Calculate how much energy escapes through the top layer
        top_emitted = boltzman_coeff * ( (sky_emissivity * (sky_temperature ** 4)) - (cell_emissivity * (T_cell ** 4)))
        top_convection = (h_air * (T_air - T_cell))

        # Calculate how much energy escapes through the bottom layer
        bottom_emitted = boltzman_coeff * ((ground_emissivity * (ground_temperature ** 4)) - (back_plate_emissivity * (T_bottom ** 4)))  # Assuming T_cell == T_bottom
        bottom_convection = (h_air * (T_air - T_bottom))  # Assuming T_cell == T_bottom

        # If the metling front approaches the end, use this if/else statement to prvent problems with very small distances
        if PCM_thickness >= x_front > upper_boundary:

            # Calculate the conduction through the PCM layer
            cell_to_alu_conduction = - shape_factor * (PCM_liq_conductivity / x_front) * (T_bottom - T_cell)

            # Calculate the gradient of the temperature of the cell
            dT_cell = panel_area * (1 / (cell_mass * cell_heat_cap)) * (cell_absorbed -
                                                                        cell_to_alu_conduction +
                                                                        top_emitted +
                                                                        top_convection)

            # Calculate the gradient of the temperature in the bottom aluminium panel
            dT_bottom = panel_area * (1 / ((PCM_liq_density * PCM_thickness * PCM_liq_heat_cap))) * (PCM_to_alu_conduction +
                                                                                                     bottom_convection +
                                                                                                     bottom_emitted)

            # Calculate the new temperatures of the cell, bottom and the glass top
            T_bottom = T_bottom + dT_bottom * dt

            # Calculate the new position and temperature of the melting front.
            x_front = PCM_thickness
            T_front = T_bottom

            # print("second: ", n + dt_ints)

            all_values[n+dt_ints] = (T_cell + dT_cell * dt, efficiency(T_cell=T_cell))

        elif lower_boundary <= x_front <= upper_boundary:

            # Calculate the conduction through the PCM layer
            cell_to_PCM_conduction = - shape_factor * (PCM_liq_conductivity / x_front) * (T_front - T_cell)
            PCM_to_alu_conduction = - PCM_sol_conductivity * (T_bottom - T_front) / (PCM_thickness - x_front)

            # Calculate the gradient of the temperature of the cell
            dT_cell = panel_area * (1 / (cell_mass * cell_heat_cap)) * (cell_absorbed -
                                                                        cell_to_PCM_conduction +
                                                                        top_emitted +
                                                                        top_convection)

            # Calculate the gradient of the temperature in the bottom aluminium panel
            dT_bottom = panel_area * (1 / ( (PCM_sol_density * (PCM_thickness - x_front) * PCM_sol_heat_cap) ) ) * (PCM_to_alu_conduction +
                                                                                                                       bottom_convection +
                                                                                                                       bottom_emitted)

            # Calculate the gradient of the temperature at the melting front
            dT_front = panel_area * (1 / (PCM_liq_density * x_front * PCM_liq_heat_cap)) * (cell_to_PCM_conduction -
                                                                                               PCM_to_alu_conduction)

            # Calculate the velocity of the melting front
            dx_front = (cell_to_PCM_conduction - PCM_to_alu_conduction) / (PCM_sol_density * PCM_latent_heat)

            # Calculate the new temperatures of the cell, bottom and the glass top
            T_bottom = T_bottom + dT_bottom * dt

            # Calculate the new position and temperature of the melting front.
            x_front = x_front + dx_front * dt
            T_front = min(T_melt, T_front + dT_front * dt)

            all_values[n + dt_ints] = (T_cell + dT_cell * dt, efficiency(T_cell=T_cell))

        elif x_front < lower_boundary:
            raise ValueError("Melting front moved backwards")

        else:
            raise ValueError("Something went wrong. Melting front position: {}m".format(x_front))

    t, values = list(zip(*all_values.items()))
    T_cell, cell_efficiency = list(zip(*values))

    mean_power = total_power / len(t)
    mean_efficiency = np.mean(np.array(cell_efficiency))
    mean_cell_temperature = np.mean(np.array(T_cell))

    if plotting == True:
        T_cell = np.array(T_cell)
        cell_efficiency = np.array(cell_efficiency)
        intervals = slice(0, len(t) + 1, int(900 / dt))
        print("Plotting...")
        plot(t=t_real[intervals], x1=T_cell[intervals], x2=(cell_efficiency[intervals]) * 100)

    return mean_power, mean_efficiency, mean_cell_temperature


if __name__ == "__main__":

    # PCM = [0.035,1490,175000.0,1.0,301.15]

    PCM_thickness, PCM_density, PCM_latent_heat, PCM_conductivity, T_melt = 0.035,1490,175000.0,1.0,301.15

    mean_power, mean_efficiency, mean_cell_temperature = simulation(plotting=True,
                                                                    PCM_thickness=PCM_thickness,
                                                                    PCM_density=PCM_density,
                                                                    PCM_latent_heat=PCM_latent_heat,
                                                                    PCM_conductivity=PCM_conductivity,
                                                                    T_melt=T_melt)

    # thicknesses = [0.045, 0.055, 0.065, 0.075]
    #
    # PCMs = [[1562, 190.8, 0.54, 302.15], [1126, 127.2, 0.189, 295.15], [878, 152.7, 0.153, 305.15], [1380, 112.0, 0.6, 294.15], [1380, 151.3, 0.6, 297.15], [1420, 162.3, 0.6, 301.15], [1420, 162.3, 0.6, 305.15], [1587, 210.0, 0.45, 319.15], [1584, 100.0, 0.43, 317.15], [2100, 115.0, 0.52, 307.15], [1460, 200.0, 0.51, 305.15], [1304, 190.0, 0.48, 303.15], [1530, 183.0, 0.54, 300.15], [1530, 180.0, 0.54, 298.15], [1530, 175.0, 0.54, 296.15], [1530, 170.0, 0.54, 295.15], [1520, 160.0, 0.43, 292.15], [1525, 160.0, 0.43, 290.15], [910, 155.0, 0.22, 319.15], [805, 242.0, 0.18, 317.15], [780, 165.0, 0.18, 316.15], [905, 105.0, 0.21, 315.15], [810, 230.0, 0.18, 313.15], [900, 105.0, 0.22, 312.15], [810, 235.0, 0.18, 310.15], [790, 217.0, 0.18, 309.15], [845, 130.0, 0.21, 305.15], [810, 225.0, 0.18, 302.15], [789, 155.0, 0.21, 301.15], [790, 150.0, 0.21, 299.15], [810, 226.0, 0.18, 298.15], [785, 150.0, 0.18, 298.15], [790, 145.0, 0.18, 297.15], [785, 145.0, 0.18, 296.15], [820, 216.0, 0.18, 295.15], [785, 145.0, 0.18, 295.15], [785, 150.0, 0.18, 290.15], [1490, 175.0, 1.0, 301.15], [1490, 175.0, 1.0, 297.15], [1490, 175.0, 1.0, 292.15], [1490, 175.0, 1.0, 290.15]]
    #
    # for PCM_thickness in thicknesses:
    #     for PCM in PCMs:
    #         PCM_density, PCM_latent_heat, PCM_conductivity, T_melt = PCM
    #         PCM_latent_heat *= 1000
    #
    #         print("Calculating for: {}, {}, {}, {}, {}".format(PCM_thickness, PCM_density, PCM_latent_heat, PCM_conductivity, T_melt))
    #
    #         mean_power, mean_efficiency, mean_cell_temperature = simulation(plotting=False,
    #                                                                         PCM_thickness=PCM_thickness,
    #                                                                         PCM_density=PCM_density,
    #                                                                         PCM_latent_heat=PCM_latent_heat,
    #                                                                         PCM_conductivity=PCM_conductivity,
    #                                                                         T_melt=T_melt)
    #
    #         if os.path.exists(os.getcwd() + "/data2.csv"):
    #             df = pandas.read_csv("data2.csv")
    #             df2 = pandas.DataFrame(data=[[PCM_thickness, PCM_density, PCM_latent_heat, PCM_conductivity, T_melt, mean_power, mean_efficiency, mean_cell_temperature]],
    #                                   columns=["PCM_thickness", "PCM_density", "PCM_latent_heat", "PCM_conductivity", "T_melt", "mean_power", "mean_efficiency", "mean_cell_temperature"])
    #             df = df.append(df2)
    #             df = df.drop_duplicates(subset=(["PCM_thickness", "PCM_density", "PCM_latent_heat", "PCM_conductivity", "T_melt"]), keep="last")
    #             df.to_csv("data2.csv", index=None)
    #         else:
    #             df = pandas.DataFrame(data=[[PCM_thickness, PCM_density, PCM_latent_heat, PCM_conductivity, T_melt, mean_power, mean_efficiency, mean_cell_temperature]],
    #                                   columns=["PCM_thickness", "PCM_density", "PCM_latent_heat", "PCM_conductivity", "T_melt", "mean_power", "mean_efficiency", "mean_cell_temperature"])
    #             df.to_csv("data2.csv", index=None)


