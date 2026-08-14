from scipy.constants import R
from scipy.integrate import odeint
import numpy as np
import pandas as pd


def A1(con):
    ca0 = 0.167  # M
    cb0 = 0.250  # M
    ar = 3.1e7  # L^(1/2) mol^(−3/2) s^(−1)
    ear = 55  # kJ/mol
    eas1 = 100  # kJ/mol
    as1 = 1e12  # s^(-1)#case 3
    eas2 = 50  # kJ/mol
    as2 = 3.1e5  # L^(1/2) mol^(−3/2) s^(−1)#case 4
    eai = {'a': 0.7, 'b': 0.4, 'c': 0.3, 'd': 0.7, 'e': 0.0, 'f': 2.2, 'g': 3.8, 'h': 7.3}  # kJ/mol
    T = con[0]
    tre = con[1]
    ccat = con[2]
    # change units
    T = T + 273.15
    tre = tre * 60
    ccat = ccat / 1000
    ea = eai[con[3]]
    kr = ccat ** 0.5 * ar * np.exp(-(ear + ea) / (T * R / 1000))

    ks2 = as2 * np.exp(-eas2 / (T * R / 1000))  # case 1

    # /1000:R J to kJ; kr:mol^(-1)s^(-1)
    def reaction(w, time):
        a, b, c, d = w
        f1 = -kr * a * b
        f2 = -kr * a * b - ks2 * b * c
        f3 = kr * a * b - ks2 * b * c
        f4 = ks2 * b * c  # case 1
        return [f1, f2, f3, f4]

    tre = tre / 10
    time = np.arange(0, tre, 0.001)
    re = odeint(reaction, (ca0, cb0, 0.0, 0.0), time)
    cr = re[-1, :][2]
    y = cr / ca0  # Reaction product yield R
    return round(y, 4)


def A2(con):
    ca0 = 0.167  # M
    cb0 = 0.250  # M
    ar = 3.1e7  # L^(1/2) mol^(−3/2) s^(−1)
    ear = 55  # kJ/mol
    eas1 = 100  # kJ/mol
    as1 = 1e12  # s^(-1)#case 3
    eas2 = 50  # kJ/mol
    as2 = 3.1e5  # L^(1/2) mol^(−3/2) s^(−1)#case 4
    eai = {'a': 0.7, 'b': 0.4, 'c': 0.3, 'd': 0.7, 'e': 0.0, 'f': 2.2, 'g': 3.8, 'h': 7.3}  # kJ/mol
    T = con[0]
    tre = con[1]
    ccat = con[2]
    # change units
    T = T + 273.15
    tre = tre * 60
    ccat = ccat / 1000
    ea = eai[con[3]]
    kr = ccat ** 0.5 * ar * np.exp(-(ear + ea) / (T * R / 1000))

    ks1 = as1 * np.exp(-eas1 / (T * R / 1000))  # case 2

    # /1000:R J to kJ; kr:mol^(-1)s^(-1)
    def reaction(w, time):
        a, b, c, d = w
        f1 = -kr * a * b
        f2 = -kr * a * b - ks1 * b
        f3 = kr * a * b
        f4 = ks1 * b  # case 2
        return [f1, f2, f3, f4]

    tre = tre / 10
    time = np.arange(0, tre, 0.001)
    re = odeint(reaction, (ca0, cb0, 0.0, 0.0), time)
    cr = re[-1, :][2]
    y = cr / ca0  # Reaction product yield R
    return round(y, 4)


def C(con):
    data = pd.read_csv('data_Pfizer.csv')
    y = data.loc[(data["Reactant_1"].str.contains(con[0], regex=False, na=False)) & (data["Reactant_2"].str.contains(con[1], regex=False, na=False)) & (
        data["Ligand"].str.contains(con[2], regex=False, na=False)) & (data["Base"].str.contains(con[3], regex=False, na=False)) & (
                           data["Solvent"].str.contains(con[4], regex=False, na=False))]['Yield']
    y = float(y.iloc[0])
    return round(y, 4)


def run_exp(data):
    reaction = data['reaction']
    n_candidates = len(data['condition']) - len(data['outcomes'])
    q = data['q']
    candidates = data['condition'][len(data['condition'])-n_candidates:][0:q]
    if 'A1' in reaction:
        for i in candidates:
            data['outcomes'].append(A1(i))
        return data
    if reaction == 'A2':
        for i in candidates:
            data['outcomes'].append(A2(i))
        return data
    if 'C' in reaction:
        for i in candidates:
            data['outcomes'].append(C(i))
        return data

