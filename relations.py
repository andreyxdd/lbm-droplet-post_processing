import math

class Relations:
  
  def __init__(self, Ca, sigma_0, R_0, Re, rho, C_sh):
    self.Ca = Ca
    self.sigma_0 = sigma_0
    self.R_0 = R_0
    self.Re = Re
    self.rho = rho
    self.C_sh = C_sh

    self.H = 8 * R_0

  def set_zeta(self, zeta):
    self.zeta = zeta

  def set_tau_g(self, tau_g):
    self.tau_g = tau_g
    
  def compute_and_set_mobility(self):
    M = self.tau_g/3
    self.M = M

  """
  From the definition of the capilary number
  """
  def compute_mu_mult_G(self):
    return self.Ca * self.sigma_0 / self.R_0

  """
  From the definition of the Reynolds number
  """
  def comput_mu_divided_by_G(self):
    return self.Re / self.rho / (self.R_0**2)

  def compute_and_set_G(self):
    A = self.compute_mu_mult_G()
    B = self.comput_mu_divided_by_G()
    self.G = math.sqrt(A * B)

    return self.G

  def compute_and_set_mu_c(self):
    A = self.compute_mu_mult_G()
    G = self.compute_and_set_G()
    self.mu_c = A / G
    
    return self.mu_c

  def compute_analytical_taylor_deformation(self):
    const_A = (35/32) * self.Ca
    const_B = (3.5/2) * self.C_sh
    const_C = (self.R_0 / self.H)**3
      
    return const_A * (1 + const_B * const_C)

  def compute_chan_number(self):
    if not self.zeta:
      raise Exception("Zeta is undefined")

    return self.zeta / self.R_0

  def compute_peclet_number(self):
    if not self.zeta:
      raise Exception("Zeta is undefined")

    if not self.M:
      raise Exception("Mobility is undefined")

    if not self.G:
      raise Exception("G is undefined")
    
    return self.G * self.R_0 * self.zeta / self.M
  
  def compute_Ch_and_Pe_from_zeta_and_tau_g(self, zeta, tau_g):
    self.set_zeta(zeta)
    self.set_tau_g(tau_g)
    self.compute_and_set_mobility()
    self.compute_and_set_G()
    self.compute_and_set_mu_c()
    
    Ch = self.compute_chan_number()
    Pe = self.compute_peclet_number()
    
    return Ch, Pe
