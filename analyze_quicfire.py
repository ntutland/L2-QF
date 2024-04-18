from pathlib import Path
from quicfire_tools import SimulationOutputs
from quicfire_tools.inputs import QUIC_fire
from matplotlib import pyplot as plt
import numpy as np


def plot_array(x, title):
    plt.figure(2)
    plt.set_cmap("viridis")
    plt.imshow(x, origin="lower")
    plt.colorbar()
    plt.title(title, fontsize=18)
    plt.show()


def get_mass_burnt(sim: SimulationOutputs, arrpath: Path, plot: bool = True):
    mburnt = sim.get_output("mburnt_integ")
    mburnt_arr = mburnt.to_numpy()
    mburnt_total = mburnt_arr[-1, 0, :, :]
    if plot:
        plot_array(mburnt_total, "percent mass burnt")
    np.savetxt(arrpath / "mass_burnt_pct.txt", mburnt_total)


def get_surface_consumption(sim: SimulationOutputs, arrpath: Path, plot: bool = True):
    dens = sim.get_output("fuels-dens")
    dens_init = dens.to_numpy(timestep=0)
    dens_final = dens.to_numpy(len(dens.times) - 1)
    surface_consumption = dens_init[0, 0, :, :] - dens_final[0, 0, :, :]
    surf_cons_pct = (
        (dens_init[0, 0, :, :] - dens_final[0, 0, :, :]) / dens_init[0, 0, :, :]
    ) * 100
    if plot:
        plot_array(surface_consumption, "surface fuel consumption")
        plot_array(surf_cons_pct, "surface consumption percent")
    np.savetxt(arrpath / "surface_consumption.txt", surface_consumption)
    np.savetxt(arrpath / "surface_consumption_pct.txt", surf_cons_pct)


def get_canopy_consumption(sim: SimulationOutputs, arrpath: Path, plot: bool = True):
    dens = sim.get_output("fuels-dens")
    dens_init = dens.to_numpy(timestep=0)
    dens_final = dens.to_numpy(len(dens.times) - 1)
    canopy_consumption = np.sum(dens_init[0, 1:, :, :], axis=0) - np.sum(
        dens_final[0, 1:, :, :], axis=0
    )
    canopy_cons_pct = (
        (
            np.sum(dens_init[0, 1:, :, :], axis=0)
            - np.sum(dens_final[0, 1:, :, :], axis=0)
        )
        / np.sum(dens_init[0, 1:, :, :], axis=0)
    ) * 100
    if plot:
        plot_array(canopy_consumption, "canopy fuel consumption")
        plot_array(canopy_cons_pct, "canopy consumption percent")
    np.savetxt(arrpath / "canopy_consumption.txt", canopy_consumption)
    np.savetxt(arrpath / "canopy_consumption_pct.txt", canopy_cons_pct)


def get_max_power(sim: SimulationOutputs, arrpath: Path, plot: bool = True):
    energy = sim.get_output("surfEnergy")
    initial = energy.to_numpy(timestep=0)
    prev = initial[0, 0, :, :]
    for t in range(1, len(energy.times)):
        temp = energy.to_numpy(timestep=t)
        temp = temp[0, 0, :, :]
        prev = np.maximum(temp, prev)
    max_power = prev.copy()
    if plot:
        plot_array(max_power, "max power")
    np.savetxt(arrpath / "max_power.txt", max_power)


def get_residence_time_from_power(
    sim: SimulationOutputs, arrpath: Path, plot: bool = True
):
    energy = sim.get_output("surfEnergy")
    initial = energy.to_numpy(timestep=0)
    prev = initial[0, 0, :, :]
    for t in range(1, len(energy.times)):
        temp = energy.to_numpy(timestep=t)
        temp = temp[0, 0, :, :]
        temp[np.where(temp > 0)] = 1
        prev = np.add(temp, prev)
    residence_time = prev.copy()
    if plot:
        plot_array(residence_time, "residence time from power")
    np.savetxt(arrpath / "residence_time_power.txt", residence_time)


def get_residence_time_from_consumption(
    sim: SimulationOutputs, arrpath: Path, plot: bool = True
):
    dens = sim.get_output("fuels-dens")
    add_to = np.zeros((sim.ny, sim.nx))
    for t in range(1, len(dens.times)):
        temp = dens.to_numpy(timestep=t)
        temp = np.sum(temp, axis=1)
        temp = temp[0, :, :]
        prev = dens.to_numpy(timestep=t - 1)
        prev = np.sum(prev, axis=1)
        prev = prev[0, :, :]
        fire_lox = np.where(temp != prev)
        no_fire_lox = np.where(temp == prev)
        temp[fire_lox] = 1
        temp[no_fire_lox] = 0
        add_to = np.add(temp, add_to)
    residence_time = add_to * 30
    if plot:
        plot_array(residence_time, "residence time from consumption")
    np.savetxt(arrpath / "residence_time_consumption.txt", residence_time)


def get_max_reaction_rate(sim: SimulationOutputs, arrpath: Path, plot: bool = True):
    react = sim.get_output("fire-reaction_rate")
    initial = react.to_numpy(timestep=0)
    initial = np.sum(initial, axis=1)
    prev = initial[0, :, :]
    for t in range(1, len(react.times)):
        temp = react.to_numpy(timestep=t)
        temp = np.sum(temp, axis=1)
        temp = temp[0, :, :]
        prev = np.maximum(temp, prev)
    max_react = prev.copy()
    if plot:
        plot_array(max_react, "max reaction rate")
    np.savetxt(arrpath / "max_reaction_rate.txt", max_react)


HERE = Path(__file__).parent
runs_dir = HERE / "7.QUICFIRE-MODEL" / "projects"

runs = ["Fire2-Dry", "Fire5-Dry", "Fire2-Wet", "Fire5-Wet"]
for run in runs:
    runpath = runs_dir / run
    quic_fire = QUIC_fire.from_file(runpath, version="v5")
    nz = quic_fire.nz
    sim_outputs = SimulationOutputs(runpath / "Output", nz=nz, ny=400, nx=400)
    print(sim_outputs.list_available_outputs())

    arrpath = runpath / "Arrays"
    arrpath.mkdir(exist_ok=True)

    print("\t- getting mass burnt")
    get_mass_burnt(sim_outputs, arrpath, False)
    print("\t- getting surface consumption")
    get_surface_consumption(sim_outputs, arrpath, False)
    print("\t- getting canopy consumption")
    get_canopy_consumption(sim_outputs, arrpath, False)
    print("\t- getting max power")
    get_max_power(sim_outputs, arrpath, False)
    print("\t- getting residence time from power")
    get_residence_time_from_power(sim_outputs, arrpath, False)
    print("\t- getting residence time from consumption")
    get_residence_time_from_consumption(sim_outputs, arrpath, False)
    print("\t- getting max reaction rate")
    get_max_reaction_rate(sim_outputs, arrpath, False)

# sim = SimulationOutputs(runpath / "Output", nz=40, ny=600, nx=1100)
# # print(sim.list_available_outputs())
# dens = sim.get_output("fuels-dens")
# dens_init = dens.to_numpy(timestep=0)
# dens_init = dens_init[0, :, :, :]
# dens_final = dens.to_numpy(timestep=len(dens.times) - 1)
# dens_final = dens_final[0, :, :, :]


# # What is the canopy base height across the domain?
# cbh_arr = np.zeros((dens_init.shape[1], dens_init.shape[2]))
# for x in range(dens_init.shape[2]):
#     for y in range(dens_init.shape[1]):
#         i = 1
#         for z in range(1, dens_init.shape[0]):
#             if dens_init[z, y, x] > 0:
#                 cbh_arr[y, x] = i
#                 break
#             else:
#                 i += 1

# plot_array(cbh_arr, "canopy base height")
