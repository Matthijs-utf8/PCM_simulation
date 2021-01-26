
# General constants
boltzman_coeff = 5.67 * (10**-8) # [W/(m^2*K^4)]
panel_area = 1 # [m^2]

# Ambience parameters
T_air = 298 # [K]
wind_speed = 10 # [m/s]
h_air = 8.91 + 2 * wind_speed # heat transfer coefficient of air on top panel [-]
solar_irridance = 1000 # [W/m^2]
sky_temperature = 0.0552 * (T_air ** 1.5)
sky_emissivity = 0.95 # [-]
ground_emissivity = 0.9 # [-]
ground_temperature = T_air

# PV cell parameters
eta_ref = 0.156 # [-]
beta_ref = 0.0045 # [-]
cell_absorptance = 0.9 # [-]
cell_thickness = 0.000225 # [m]
cell_density = 2330 # [kg/m^3]
cell_heat_cap = 677 # [J/(kg*K)]
cell_emissivity = 0.9 # [-]
cell_mass = cell_density * cell_thickness * panel_area # [kg]

# Back plate parameters
back_plate_emissivity = 0.9 # [-]

# Time parameters
t_start = 0 # [s]
t_end = 8*60*60 # [s]
dt = 0.01 # [s]
factor = int(1/dt)
t_end_ints = t_end * factor
dt_ints=factor*dt



