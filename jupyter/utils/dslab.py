# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:47:39 2021

@author: Adolfo
"""

import numpy as np

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

def Vb(b, sigma, pc, ps, m):
  return 1/(2*np.sqrt(1-b))*(m*np.pi + np.arctan(pc*np.sqrt(b+sigma)/np.sqrt(1-b)) + np.arctan(ps*np.sqrt(b)/np.sqrt(1-b)))  

def b2neff(ns, NA, b):
  return np.sqrt(b*NA**2+ns**2)