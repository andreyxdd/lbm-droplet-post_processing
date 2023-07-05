import os, shutil
from pprint import pprint

import numpy as np
import relations
import matplotlib.pyplot as plt

# --- Settings
round_to = 6 # decimals
skip_timseries_plot = False
plot_timeseries = False
plot_contour = False

zeta_array = [2, 3, 5, 6]
mobilities_array = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] #, 2, 3]
tau_g_array = [float(round(3*mobility, 3)) for mobility in mobilities_array]
# ---

if os.path.exists('./plots_dt') and os.path.isdir('./plots_dt'):
  shutil.rmtree('./plots_dt')
os.mkdir('plots_dt')

if __name__ == '__main__':
  r = relations.Relations(
    Ca = 0.1,
    sigma_0 = 0.01,
    R_0 = 30,
    Re = 0.0625,
    rho = 1,
    C_sh = 5.6996
  )

  D_t = r.compute_analytical_taylor_deformation()  
  print(f"Taylor deformation (analytical solution) = {D_t}")

  heat_map = {}
  for idx_zeta, zeta in enumerate(zeta_array):
    for idx_tau, tau_g in enumerate(tau_g_array):
      Ch, Pe = r.compute_Ch_and_Pe_from_zeta_and_tau_g(zeta, tau_g)
      Ch_rouned = round(Ch, round_to)
      Pe_rounded = round(Pe, round_to)

      results_dt = []
      with open(f'./results_dt/zeta-{zeta}__tau_g-{tau_g}') as f:
        for line in f:
            results_dt.append(float(line))

      M = mobilities_array[idx_tau]

      # maximum D_T
      max_D_T = max(results_dt)
      
      # filling heat map
      heat_map[(idx_zeta, idx_tau)] = (Ch_rouned, Pe_rounded, round(max_D_T, round_to))

      # plot D_T over time for the given zeta and tau_g
      if not skip_timseries_plot:
        fig, ax = plt.subplots()
        ax.plot(results_dt, label =  r'Numerical')
        plt.title(
          f"Ch = {Ch_rouned}, Pe = {Pe_rounded}"
          +
          f" (zeta = {zeta}, M = {M})"
        )
        plt.xlabel('Time')
        plt.ylabel(r'Taylor deformation, $D_T$')
        
        # numerical D_T
        plt.axhline(y = max_D_T, linestyle = '--')
        
        # analytical D_T
        plt.axhline(y = D_t, linestyle = '-.', label =  'Analytical', color = 'tab:orange')  

        plt.legend()

        if plot_timeseries:
          plt.show()
        else:
          fig.savefig(f'./plots_dt/zeta-{zeta}__mobility-{M}.png')

  # --- post-process
  x = []
  y = []
  z = np.zeros((len(zeta_array), len(mobilities_array)))
  for ((i,j), v) in heat_map.items():
    x.append(v[0])
    y.append(v[1])
    z[i, j] = v[2]

  x_max = max(x)
  x_min = min(x)
  y_max = max(y)
  y_min = min(y)
  # ---

  # --- prepare for plotting
  # X, Y = np.meshgrid(np.array(list(set(x))), np.array(list(set(y))))
  X, Y = np.meshgrid(
    np.linspace(x_min, x_max, len(zeta_array)),
    np.linspace(y_min, y_max, len(mobilities_array)),
  )
  Z = z.transpose()
  
  # pprint(X)
  # pprint(Y)
  # pprint(Z)
  # ---

  # --- plot contour
  plt.close('all')
  
  fig, ax = plt.subplots(1, 1)
  C = ax.contourf(X, Y, Z)
  ax.set_title(r'Taylor deformation, $D_T$')
  ax.set_xlabel('Ch')
  ax.set_ylabel('Pe')
  fig.colorbar(C)

  if plot_contour:
    plt.show()
  else:
    fig.savefig('D_t-contour')
  # ---
