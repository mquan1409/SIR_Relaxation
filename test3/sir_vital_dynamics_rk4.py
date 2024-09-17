import math
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

def g_hat(r: float, p: int, args: dict) -> float:
    if r < 0:
        print('--r:',r)
    muy = args['beta'] / args['gamma']
    delta_t = args['T'] / args['P']
    result = args['gamma']*args['n']*math.exp(muy*(math.exp((-1)*args['sigma']*p*delta_t) - 1))\
            *math.exp( (-1)*muy*math.exp( (-1)*args['sigma']*p*delta_t )*r ) \
            + args['gamma']*r
    return result

def N(p:int, args) -> float:
    delta_t = args['T'] / args['P']
    result = math.exp((-1)*args['sigma']*p*delta_t)*args['N0']
    return result

def K1(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_hat']*R[k][p-1] + args['gamma']*args['N0'] - g_hat(R[k-1][p-1],(p - 1),args) + args['M_hat']*R[k-1][p-1])
    return result

def K2(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_hat']*(R[k][p-1] + K1(tp,p,k,R,args)/2) + args['gamma']*args['N0']\
            - g_hat((1/2)*(R[k-1][p-1] + R[k-1][p]),(p - 1/2),args) + args['M_hat']*(1/2)*(R[k-1][p-1] + R[k-1][p]))
    return result

def K3(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_hat']*(R[k][p-1] + K2(tp,p,k,R,args)/2) + args['gamma']*args['N0']\
            - g_hat((1/2)*(R[k-1][p-1] + R[k-1][p]),(p - 1/2),args) + args['M_hat']*(1/2)*(R[k-1][p-1] + R[k-1][p]))
    return result

def K4(tp: float, p: int, k: int, R, args: dict):
    delta_t = args['T'] / args['P']
    result = delta_t * ((-1)*args['M_hat']*(R[k][p-1] + K3(tp,p,k,R,args)) + args['gamma']*args['N0']\
            - g_hat(R[k-1][p],p,args) + args['M_hat']*R[k-1][p])
    return result

def Rpk_bar(args: dict,R_values,S_values,I_values,days) -> float:
    TOL = 0.01
    stop = False
    R_bar = np.zeros((args['K'] + 1,args['P'] + 1))
    delta_t = args['T'] / args['P']
    I_value_max = 0
    I_value_max_day = 0
    for k in range(1,args['K'] + 1):
        for p in range(1,args['P'] + 1):
            R_bar[k][p] = R_bar[k][p-1] + (1/6)*K1(p*delta_t,p,k,R_bar,args) + (1/3)*K2(p*delta_t,p,k,R_bar,args)\
                    + (1/3)*K3(p*delta_t,p,k,R_bar,args) + (1/6)*K4(p*delta_t,p,k,R_bar,args)
            if k == args['K']:
                R_value = R_bar[k][p]/math.exp(args['sigma']*p*delta_t)
                S_value = S(R_value,p,args)
                I_value = I(R_value,S_value,p,args)
                R_values.append(R_value)
                S_values.append(S_value)
                I_values.append(I_value)
                days.append(p*delta_t)
                if (I_value > I_value_max):
                    I_value_max = I_value
                    I_value_max_day = p*delta_t            
    
    for k in range(2,args['K'] + 1):
        if stop:
            break
        for p in range(2,args['P'] + 1):
            R_old = R_bar[k-1][p-1]/math.exp(args['sigma']*p*delta_t)
            R_new = R_bar[k][p]/math.exp(args['sigma']*p*delta_t)
            diff = abs(R_new - R_old)/abs(R_old)
            if diff < TOL:
                print('Satisfied K P')
                print('--- K:',k,',P:',p)
                stop = True 
                break

    print('--- I value max =', I_value_max, 'on day', I_value_max_day)

    return R_bar[args['K']][args['P']]

def S(R: float,p,args: dict) -> float:
    delta_t = args['T'] / args['P']
    muy = args['beta'] / args['gamma']
    S_bar = args['n'] * math.exp( muy * ( math.exp((-1)*args['sigma']*p*delta_t) - 1))\
            * math.exp((-1)*muy*R)
    result = S_bar / math.exp(args['sigma']*p*delta_t)
    return result

def I(R: float, S: float, p: int, args: dict) -> float:
    result = N(p,args) - R - S
    return result

def main() -> None:
    R_values_old = None
    args = {
            'N0':1000,
            'M_hat':0.02,
            'sigma':0.001,
            'gamma':0.02,
            'beta':0.0004,
            'T':365,
            'P':100,
            'K':5,
            'n':998
            }
    R_values = []
    S_values = []
    I_values = []
    days = []
    print('--- K:',args['K'],',P:',args['P'])
    Rpk_bar(args,R_values,S_values,I_values,days)
    '''
    muy = args['beta']/args['gamma']
    S_when_I_max = ((args['gamma']+args['sigma'])/args['beta'])
    I_value_max_true = S_when_I_max*(math.log(S_when_I_max) + args['sigma']*args['T'] - 1 - math.log(n))\
            + (args[N0] - args[n]) + args[n]
    if args['alpha'] >= args['sigma']:
        I_value_max_true += args['N0']*(math.exp((args['alpha'] - args['sigma'])*args['T']) - math.exp((-1)*args['sigma']*args['T']))
    if args['alpha'] < args['sigma']:
        I_value_max_true += args['N0']*(1 - math.exp((-1)*muy*args['T']))
    '''
    plt.figure(figsize=(20,12))
    plt.plot(days,R_values,label='R',linewidth=4.5)
    plt.plot(days,S_values,label='S',linewidth=4.5)
    plt.plot(days,I_values,label='I',linewidth=4.5)
    plt.title('SIR model ('\
            + 'gamma: ' + str(args['gamma'])\
            + ', beta: ' + str(args['beta'])\
            + ', sigma: ' + str(args['sigma'])\
            + ', M_hat: ' + str(args['M_hat'])\
            + ', P: ' + str(args['P'])\
            + ', K: ' + str(args['K']) + ')', fontsize=25)
    plt.xlabel('Days', fontsize=25)
    plt.ylabel('People', fontsize=25)
    plt.tick_params(labelsize=25)
    plt.legend(prop={ 'size': 40 }, loc='center right', fontsize=25)
    plt.savefig('sir_vital_dynamics_rk4.png')
    print('--- plotted ' + 'P: ' + str(args['P']) + ', K: ' + str(args['K']) + ' ---')
    plt.show()
    plt.clf()


if __name__=='__main__':
    #start_time = time.time()
    main()
    #stop_time = time.time()
    #execution_time = stop_time - start_time
    #print('Execution time: ',execution_time)
