import numpy as np
import matplotlib.pyplot as plt

# Draw masses from the maschberger IMF as defined by the equation 4 in
# table 1 of the paper Maschberger T., 2013, MNRAS, 429, 1725 "On the 
# function describing the stellar initial mass function"
def maschberger_IMF(N_stars,  m_lower_lim, m_upper_lim, alpha=2.3, beta=1.4, mu=0.2):

    # Set up shorthand 
    alpha_minus = 1. - alpha
    beta_minus = 1. - beta
    
    # Get G as a function of the desired upper and lower mass linits
    G_m_lower_lim = (1. + (m_lower_lim / mu)**(alpha_minus))**(beta_minus)
    G_m_upper_lim = (1. + (m_upper_lim / mu)**(alpha_minus))**(beta_minus)
    
    # Draw a number from a uniform distribution for each star
    u = np.random.uniform(size=N_stars)
    
    # Get G(m) for each random draw
    G_m = u * (G_m_upper_lim - G_m_lower_lim) + G_m_lower_lim
    
    # Convert G(m) into an array holding the stellar masses
    m = mu*(G_m**(1./beta_minus) - 1.)**(1./alpha_minus)
    
    return m


def plot_loghist(x, bins=30):
  hist, bins = np.histogram(x, bins=bins)
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  plt.hist(x, bins=logbins)
  plt.xscale('log')
  plt.yscale('log')
  return

def plot_IMF_dist(N_stars=100000, m_lower_lim=0.01, m_upper_lim=100):

    m = maschberger_IMF(N_stars, m_lower_lim, m_upper_lim)
    
    plot_loghist(m, bins=30)
    plt.show()
    
    return

