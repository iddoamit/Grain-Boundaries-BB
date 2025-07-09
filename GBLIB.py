#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 17:00:16 2025

@author: cmrv41
"""

import numpy as np
import pandas as pd
from scipy.optimize import root
import subprocess
import ltspice
import os

here = os.path.dirname(os.path.realpath(__file__))+'\\'

'''
Constants
'''
q = 1.6e-19
e0 = 8.85e-14
me = 9.1e-31
meff = 0.86
ni = 1e10
k = 8.6e-5
mu = 1450

'''
Parameters
'''
params = pd.read_excel(here+'params.xlsx', sheet_name='Materials')
N = params['N'][0]*np.array([1., 1.])
er = params['er'][0]
NT = params['NT'][0]
T = params['T'][0]
A = params['A'][0]
L_nom = params['L_nom'][0]
M = 600
NC = 2.5e19*(meff**1.5)

'''
Electrostatics
'''
def x1_calc(x2, N=N, NT=NT):
    x1 = (NT-N[1]*x2)/N[0]
    return x1

def potential_balance(x2, VA, N=N, NT=NT):
    x1 = x1_calc(x2, N, NT)
    A = q*N[0]*(x1**2)/(2*e0*er)
    B = q*N[1]*(x2**2)/(2*e0*er)
    C = k*T*np.log(N[0]/N[1])
    return A - B - C + VA

def x2_calc(x2i, VA, N=N, NT=NT):
    S = root(potential_balance, x2i, args=(VA))
    return S.x[0]

def Vbi_calc(N=N, NT=NT):
    X = np.array([0., 0.])
    X[1] = x2_calc(NT/(2*np.mean(N)), 0)
    X[0] = x1_calc(X[1], N, NT)
    Vbi = q*N*(X**2)/(2*e0*er)
    return Vbi

'''
Current
'''

def define_VA(N=N, NT=NT):
    Vbi = Vbi_calc(N, NT)
    VA = np.linspace(-1.5*np.round(np.max(Vbi), 1), 1.5*np.round(np.max(Vbi), 1), M)
    return VA
    
def calc_dep(VA, N=N, NT=NT):
    X = np.zeros((M, 2))
    x2i = NT/(2*np.mean(N))
    for m in range(M//2, M):
        if m == M//2:
            X[m, 1] = x2_calc(x2i, VA[m], N, NT)
        else:
            X[m, 1] = x2_calc(X[m-1, 1], VA[m], N, NT)
        X[m, 0] = x1_calc(X[m, 1], N, NT)
    for m in range(M//2):
        X[M//2-m-1, 1] = x2_calc(X[M//2-m, 1], VA[M//2-m-1], N, NT)
        X[M//2-m-1, 0] = x1_calc(X[M//2-m-1, 1], N, NT)
    return X

    
    
def calc_barrier(X, N=N):
    ECEF = k*T*np.log(NC/N)
    Phi = q*N*(X**2)/(2*e0*er)+ECEF
    return Phi
    
def calc_current(Phi):
    TE = np.sum(np.exp(-Phi/(k*T))*np.array([1, -1]), 1)
    AA = 120*meff
    I = A*AA*(T**2)*TE
    return I
    
def calc_resistance(VA, I):
    return VA/I
    
def resistance(N=N, NT=NT):
    VA = define_VA(N, NT)
    X = calc_dep(VA, N, NT)
    Phi = calc_barrier(X, N)
    I = calc_current(Phi)
    R = calc_resistance(VA, I)
    output = (np.vstack((VA, R))).T
    return output, X

def sig_calc(N=N):
    return q*mu*np.mean(N)

'''
Netlist
'''

simpar = pd.read_excel(here+'params.xlsx', sheet_name='Simulation')

def make_Lvec(simpar):
    Lvec = L_nom*np.random.normal(1, np.sqrt(simpar['Grain_var'][0]), 3)
    return Lvec

def gen_GB_txt(output, node):
    txt = f'R=table(V(n{node:03.0f},n{node+1:03.0f}),'
    for m in range(len(output)):
        txt += f'{output[m, 0]:.3f},'
        txt += f'{output[m, 1]:.4e},'
    txt = txt[:-1]
    txt += ')'
    return txt

def gen_GL_txt(output, X, sig, node, gn, Lvec):
    qnr = Lvec[gn]-X[:, 0]
    qnr_R = qnr/(A*sig)
    txt = f'R=table(V(n{node+1:03.0f},n{node+2:03.0f}),'
    for m in range(len(output)):
        txt += f'{output[m, 0]:.3f},'
        txt += f'{qnr_R[m]:.4e},'
    txt = txt[:-1]
    txt += ')'
    return txt

def gen_GR_txt(output, X, sig, node, gn, Lvec):
    qnr = Lvec[gn]-X[:, 1]
    qnr_R = qnr/(A*sig)
    txt = f'R=table(V(n{node-1:03.0f},n{node:03.0f}),'
    for m in range(len(output)):
        txt += f'{output[m, 0]:.3f},'
        txt += f'{qnr_R[m]:.4e},'
    txt = txt[:-1]
    txt += ')'
    return txt

def write_netlist(simpar, output, X, sig):
    ff = open(here+simpar['filename'][0]+'.net', 'w')
    Lvec = make_Lvec(simpar)
    
    # Simulation title
    txt = f"* A grain size variate simulation with variation of {simpar['Grain_var'][0]*100:.0f}%, and grain sizes of "
    for i in range(3):
        txt += f"{1e4*Lvec[i]:.2f}, "
    txt = txt[:-2]
    txt += 'um'+'\n'
    ff.write(txt)
    
    # Bias source
    ff.write('V1 n001 0 1'+'\n')
    
    #Left most half grain
    R = Lvec[0]/(2*A*sig)
    ff.write(f"R1L n001 n002 {R:.4e}"+'\n')
    
    #Left dynamic half
    R = gen_GL_txt(output, X, sig, 2, 0, Lvec)
    ff.write('R1R n002 n003 '+R+'\n')
    
    #First grain boundary
    R = gen_GB_txt(output, 3)
    ff.write('RGB1 n003 n004 '+R+'\n')
    
    #Left half middle grain
    R = gen_GR_txt(output, X, sig, 4, 1, Lvec)
    ff.write('R2L n004 n005 '+R+'\n')
    
    #Right half middle grain
    R = gen_GL_txt(output, X, sig, 5, 1, Lvec)
    ff.write('R2R n005 n006 '+R+'\n')
    
    #Second grain boundary
    R = gen_GB_txt(output, 3)
    ff.write('RGB2 n006 n007 '+R+'\n')
    
    #Right dynamic half
    R = gen_GR_txt(output, X, sig, 7, 2, Lvec)
    ff.write('R3L n007 n008 '+R+'\n')
    
    #Right most half grain
    R = Lvec[2]/(2*A*sig)
    ff.write(f"R3R n008 0 {R:.4e}"+'\n')
    
    #Simulation directive
    ff.write(f'.dc V1 {3*output[0, 0]:.2f} {3*output[-1, 0]:.2f} 0.005'+'\n')
    ff.write('.backanno'+'\n')
    ff.write('.end')
    
    ff.close()

def run_simulation(simpar):
    subprocess.run([simpar['LT_Path'][0], '-b', here+simpar['filename'][0]+'.net'])

def analyse(simpar):
    LT = ltspice.Ltspice(here+simpar['filename'][0]+'.raw')
    LT.parse()
    VA = LT.get_data('V(n001)')
    I = -LT.get_data('I(V1)')
    IV = (np.vstack((VA, I))).T
    return IV

def run_batch(simpar):
    n = simpar['Iterations'][0]
    output, X = resistance(N, NT)
    sig = sig_calc(N)
    
    for i in range(n):
        print(f"Simulation iteration number {i+1:.0f} of {simpar['Iterations'][0]:.0f}")
        write_netlist(simpar, output, X, sig)
        run_simulation(simpar)
        O = analyse(simpar)
        if i == 0:
            data = np.copy(O)
        else:
            iO = np.interp(data[:, 0], O[:, 0], O[:, 1])
            data = np.column_stack((data, iO))
    return data

'''
MAIN
'''
def main():
    simpar = pd.read_excel(here+'params.xlsx', sheet_name='Simulation')
    data = run_batch(simpar)
    R = np.zeros((data.shape[0], data.shape[1]-1))
    for i in range(R.shape[1]):
        R[:, i] = data[:, 0]/data[:, i+1]
    col_I = ['Voltage (V)']
    col_R = ['Voltage (V)']
    for i in range(R.shape[1]):
        col_I.append(f'Exp, {i+1:.0f}, Curr. (A)')
        col_R.append(f'Exp, {i+1:.0f}, Res. (Ohm)')
    DATA = pd.DataFrame(data, columns=col_I)
    RES = pd.DataFrame(np.column_stack((data[:, 0], R)), columns=col_R)
    with pd.ExcelWriter(here+simpar['filename'][0]+'.xlsx') as writer:
        DATA.to_excel(writer, sheet_name='Current', index=None)
        RES.to_excel(writer, sheet_name='Resistance', index=None)
    
        
print('Starting simulation')

main()

print('Simulation finished.')
print('Results have been saved to '+here+simpar['filename'][0]+'.xlsx')
