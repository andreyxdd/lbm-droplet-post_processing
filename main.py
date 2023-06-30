import relations

zeta_array = [2, 3] # 5, 6]
mobilities_array = [0.1, 0.15, 0.2, 0.3, 0.4] #, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5]
tau_g_array = [round(mobility/3, 3) for mobility in mobilities_array]

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

  for zeta in zeta_array:
    for tau_g in tau_g_array:
      r.set_zeta(zeta)
      r.set_tau_g(tau_g)
      r.compute_and_set_mobility()
      r.compute_and_set_G()
      r.compute_and_set_mu_c()
      
      Ch = r.compute_chan_number();
      print(f"Ch number = {Ch}")

      Pe = r.compute_peclet_number()
      print(f"Peclet number = {Pe}")
