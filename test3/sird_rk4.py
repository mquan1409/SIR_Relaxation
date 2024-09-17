import math
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

def g_bar(r: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    result = args['gamma']*args['n']*math.exp((-1)*muy*r) \
           + (args['gamma'] + args['sigma'])*r
    return result

def K1(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_bar']*R[k][p-1] \
           + args['gamma']*args['N'] \
           - g_bar(R[k-1][p-1],args) + args['M_bar']*R[k-1][p-1])
    return result

def K2(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_bar']*(R[k][p-1] \
           + K1(tp,p,k,R,args)/2) \
           + args['gamma']*args['N'] \
           - g_bar((1/2)*(R[k-1][p-1] + R[k-1][p]),args) \
           + args['M_bar']*(1/2)*(R[k-1][p-1] + R[k-1][p]))
    return result

def K3(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_bar']*(R[k][p-1] \
           + K2(tp,p,k,R,args)/2) \
           + args['gamma']*args['N'] \
           - g_bar((1/2)*(R[k-1][p-1] + R[k-1][p]),args) \
           + args['M_bar']*(1/2)*(R[k-1][p-1] + R[k-1][p]))
    return result

def K4(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_bar']*(R[k][p-1] \
           + K3(tp,p,k,R,args)) + args['gamma']*args['N'] \
           - g_bar(R[k-1][p],args) + args['M_bar']*R[k-1][p])
    return result

def D(R: float, args: dict) -> float:
    result = (args['sigma']/args['gamma'])*R
    return result

def S(R: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    result = args['n'] * math.exp((-1) * muy * R)
    return result

def I(R: float, D: float, S: float, args: dict) -> float:
    result = args['N'] - R - S - D
    return result

def Rpk(args: dict,R_values,D_values,S_values,I_values,days) -> float:
    R = np.zeros((args['K'] + 1,args['P'] + 1))
    delta_t = args['T'] / args['P']
    I_value_max = 0
    I_value_max_day = 0
    for k in range(1,args['K'] + 1):
        for p in range(1,args['P'] + 1):
            R[k][p] = R[k][p-1] \
                    + (1/6)*K1(p*delta_t,p,k,R,args) \
                    + (1/3)*K2(p*delta_t,p,k,R,args) \
                    + (1/3)*K3(p*delta_t,p,k,R,args) \
                    + (1/6)*K4(p*delta_t,p,k,R,args)
            if k == args['K']:
                R_value = R[k][p]
                D_value = D(R_value,args)
                S_value = S(R_value,args)
                I_value = I(R_value,D_value,S_value,args)
                R_values.append(R_value)
                D_values.append(D_value)
                S_values.append(S_value)
                I_values.append(I_value)
                days.append(p*delta_t)
                if (I_value > I_value_max):
                    I_value_max = I_value
                    I_value_max_day = p*delta_t

    return (I_value_max,I_value_max_day)


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
    I_value_max_prev = 0.1
    I_value_max_day_prev = 0.1
    R_values = []
    D_values = []
    S_values = []
    I_values = []
    days = []

    R_values.clear()
    D_values.clear()
    S_values.clear()
    I_values.clear()
    (I_value_max,I_value_max_day) = \
        Rpk(args,R_values,D_values,S_values,I_values,days) 
    I_max_variance = \
        abs((I_value_max - I_value_max_prev)/I_value_max_prev) * 10**6
    I_max_day_variance = \
        abs((I_value_max_day - I_value_max_day_prev)/I_value_max_day_prev) * 10**3
    print('Imax: ', I_value_max, 'day', I_value_max_day)
    I_value_max_prev = I_value_max
    I_value_max_day_prev = I_value_max_day

    muy = args['beta']/args['gamma']
    intermidiate_value = ((args['gamma'] + args['sigma'])/args['beta'])
    I_value_max_true = intermidiate_value*math.log(intermidiate_value) \
                     - intermidiate_value \
                     + (args['N'] - args['n'] - args['d']) \
                     + args['n'] - intermidiate_value*math.log(args['n'])
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
    plt.savefig('sird_rk4.png')
    plt.show()
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
