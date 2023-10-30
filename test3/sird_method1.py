import math
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

def g_bar(r: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    result = args['gamma']*args['n']*math.exp((-1)*muy*r) + (args['gamma'] + args['sigma'])*r
    return result

def D(R: float, args: dict) -> float:
    result = (args['sigma']/args['gamma'])*R
    return result

def S(R: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    result = args['n'] * math.exp((-1) * muy * R)
    return result

def I(R: float, S: float, D: float, args: dict) -> float:
    result = args['N'] - R - S - D
    return result

def Rpk(args: dict,R_values,D_values,S_values,I_values,days) -> float:
    R = np.zeros((args['K'] + 1,args['P'] + 1))
    delta_t = args['T'] / args['P']
    I_value_max = 0
    I_value_max_day = 0
    for k in range(1,args['K'] + 1):
        for p in range(1,args['P'] + 1):
            R[k][p] = (R[k][p-1] + delta_t * (args['gamma']*args['N'] - g_bar(R[k-1][p],args) + args['M_bar']*R[k-1][p])) \
                    / (1 + delta_t*args['M_bar'])
            if k == args['K']:
                R_value = R[k][p]
                D_value = D(R_value,args)
                S_value = S(R_value,args)
                I_value = I(R_value,S_value,D_value,args)
                R_values.append(R_value)
                D_values.append(D_value)
                S_values.append(S_value)
                I_values.append(I_value)
                days.append(p*delta_t)
                if (I_value > I_value_max):
                    I_value_max = I_value
                    I_value_max_day = p*delta_t
    print('--- I value max =', I_value_max, 'on day', I_value_max_day)

    return R[args['K']][args['P']]


def main() -> None:
    R_values_old = None
    args = {
            'N':1000,
            'M_bar':0.015,
            'gamma':0.02,
            'beta':0.0004,
            'sigma':0.01,
            'T':365,
            'P':100,
            'K':5,
            'n':998,
            'd':0
            }
    lessThanTolerance = False
    tolerance = 0.1
    R_values = []
    D_values = []
    S_values = []
    I_values = []
    days = range(args['T'] + 1)
    days = []

    R_values.clear()
    S_values.clear()
    I_values.clear()
    D_values.clear()
    Rpk(args,R_values,D_values, S_values,I_values,days)
    intermidiate_value = ((args['gamma'] + args['sigma'])/args['beta'])
    I_value_max_true = intermidiate_value*math.log(intermidiate_value) - intermidiate_value \
            + (args['N'] - args['n'] - args['d']) + args['n'] - intermidiate_value*math.log(args['n'])
    plt.figure(figsize=(20,12))
    plt.plot(days,R_values,label='R',linewidth=4.5)
    plt.plot(days,D_values,label='D',linewidth=4.5)
    plt.plot(days,S_values,label='S',linewidth=4.5)
    plt.plot(days,I_values,label='I',linewidth=4.5)
    plt.title('SIRD model ('\
            + 'gamma: ' + str(args['gamma'])\
            + ', beta: ' + str(args['beta'])\
            + ', M_bar: ' + str(args['M_bar'])\
            + ', P: ' + str(args['P'])\
            + ', K: ' + str(args['K']) + ')', fontsize=25)
    plt.xlabel('Days', fontsize=25)
    plt.ylabel('People', fontsize=25)
    plt.tick_params(labelsize=25)
    plt.legend(prop = { "size": 40 }, loc='center right')
    plt.savefig('sird_method1.png')
    print('--- plotted ' + 'P: ' + str(args['P']) + ', K: ' + str(args['K']) + ' ---') 
    print('--- I value max true =', I_value_max_true)
    plt.show()
    plt.clf()


if __name__=='__main__':
    #start_time = time.time()
    main()
    #stop_time = time.time()
    #execution_time = stop_time - start_time
    #print('Execution time: ',execution_time)
