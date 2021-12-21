# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:47:39 2021

@author: Adolfo

Reference:
1 - Rohan D. Kekatpure, Aaron C. Hryciw, Edward S. Barnard, and Mark L. Brongersma, 
"Solving dielectric and plasmonic waveguide dispersion relations on a pocket calculator," 
Opt. Express 17, 24112-24129 (2009) """

import numpy as np
from numpy.lib import scimath as SM

def adguide(nf, nc, ns, lmbd, a, mode):
  mmax = 1000
  tol = 1e-10
  r = 0.3

  if (mode=='TM'):
    pc = (nf/nc)**2
    ps = (nf/ns)**2
  else:
    pc = 1
    ps = 1

  delta = (ns**2-nc**2)/(nf**2-ns**2)
  k0 = 2*np.pi/lmbd
  NA = np.sqrt(nf**2-ns**2)
  R = k0*a*NA

  M = int(np.ceil((2*R-np.arctan(pc*np.sqrt(delta)))/np.pi))
  m = np.arange(0,M)

  u = R*np.ones(M)
  v = np.zeros(M)
  w = R*np.sqrt(delta)*np.ones(M)

  Nit = 1

  while True:
    F = 0.5*(np.pi*m + np.arctan(ps*v/u) + np.arctan(pc*w/u))
    u_new = r*F + (1-r)*u
    if np.any(np.abs(u_new-u) <= tol):
        break
    Nit = Nit + 1
    u = u_new
    v = np.sqrt(R**2-u**2)
    w = np.sqrt(delta*R**2+v**2)
    if Nit > mmax:
        break

  kf = u/a
  gamma_s = v/a
  gamma_c = w/a
  be = np.sqrt((nf*k0)**2 - kf**2)
  return be,be/k0, kf, gamma_s, gamma_c, M

def mercatilli(nf, nc, ns, lmbd, h, mode, m, r):
  mmax = 1000
  tol = 1e-10

  if (mode=='TM'):
    pc = (nf/nc)**2
    ps = (nf/ns)**2
  else:
    pc = 1
    ps = 1

  delta = (ns**2-nc**2)/(nf**2-ns**2)
  k0 = 2*np.pi/lmbd
  NA = np.sqrt(nf**2-ns**2)
  R = k0*h*NA

  if m == 'none':
    M = int(np.ceil((2*R-np.arctan(pc*np.sqrt(delta)))/np.pi))
    m = np.arange(0,M)

  gamma_c = SM.sqrt((k0*(nf-nc))**2)
  #gamma_s = SM.sqrt((k0*(nf-ns))**2)
  #Gc = gamma_c
  #Gs = gamma_s
  k   = gamma_c
  Kc = k0*SM.sqrt(nf**2-nc**2)
  #Ks = k0*SM.sqrt(nf**2-ns**2)

  Nit = 1
  while True:
    #k_new = 2/a*(m * np.pi + np.arctan((pc*ps*gamma_c*gamma_s-k**2)-Gc*Gs)/(k*(pc*gamma_c+ps*gamma_s)))    
    #k_new = ((pc*ps*Kc*Ks)**2-(Gc*Gs)**2*(np.cos(k*h))**2+((pc*ps)**2-1)*k**4)
    #k_new = SM.sqrt(k_new/(((pc*ps)**2)*(Kc**2+Ks**2)+2*Gc*Gs*np.cos(k*h)))
    k_new = Kc/(SM.sqrt(1+(pc**-2)*(np.tan(k*h/2))**2))
    k_new = r*k_new + (1-r)*k
    if np.any(np.abs(k_new-k) <= tol):
        break
    Nit = Nit + 1
    k = k_new
    gamma_c = SM.sqrt((k0*(nf-nc))**2 - k**2)
    #gamma_s = SM.sqrt((k0*(nf-ns))**2 - k**2)
    #Gc = SM.sqrt(k**2+(pc*gamma_c)**2)
    #Gs = SM.sqrt(k**2+(ps*gamma_s)**2)
    if Nit > mmax:
        break

  return k

def Vb(b, sigma, pc, ps, m):
  return 1/(2*np.sqrt(1-b))*(m*np.pi + np.arctan(pc*np.sqrt(b+sigma)/np.sqrt(1-b)) + np.arctan(ps*np.sqrt(b)/np.sqrt(1-b)))  

def b2neff(ns, NA, b):
  return np.sqrt(b*NA**2+ns**2)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return array[idx]