import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat



def plot_flow_quantities(output_fp, r_out, u_out):
    # Save data
    rho_np = r_out.to_numpy()
    u_np = u_out.to_numpy()
    savemat(output_fp, {'r_out': rho_np, 'u_out': u_np})

    # Initialize matplotlib figure
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    cf1 = axs[0].contourf(rho_np, levels=20, cmap='viridis')
    axs[0].set_title("Density")
    fig.colorbar(cf1, ax=axs[0])

    cf2 = axs[1].contourf(u_np[:,:,0], levels=20, cmap='viridis')
    axs[1].set_title("Ux")
    fig.colorbar(cf2, ax=axs[1])

    cf3 = axs[2].contourf(u_np[:,:,1], levels=20, cmap='viridis')
    axs[2].set_title("Uy")
    fig.colorbar(cf3, ax=axs[2])

    plt.tight_layout()
    plt.show()


def generate_gif(u_mat_frames):
    savemat(".\\output.mat", {'u_mat_frames': u_mat_frames})


