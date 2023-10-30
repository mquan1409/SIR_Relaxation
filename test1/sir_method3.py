import math
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

def Ra(t: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    eta = math.sqrt(2*args['n']*muy*muy*(args['N'] - args['n']) + (args['n']*muy - 1)**2)
    psi = math.atanh((1/eta) * (args['n']*muy - 1))
    R = (1/(args['n']*muy*muy)) * (args['n']*muy - 1 + eta*math.tanh((args['gamma']*eta*t)/2 - psi))
    return R

def S(R: float, args: dict) -> float:
    muy = args['beta'] / args['gamma']
    result = args['n'] * math.exp((-1) * muy * R)
    return result

def I(R: float, S: float, args: dict) -> float:
    result = args['N'] - R - S
    return result

def main() -> None:
    args = {
            'N':1000,
            'M':0,
            'gamma':0.02,
            'beta':0.0004,
            'T':365,
            'P':1000,
            'K':0,
            'n':998
            }
    R_values = []
    S_values = []
    I_values = []
    I_value_max = 0
    I_value_max_day = 0
    days = range(args['T'] + 1)

    I_value_max = 0
    I_value_max_day = 0
    R_values.clear()
    S_values.clear()
    I_values.clear()
    for j in days:
        R_value = Ra(j,args)
        #if R_value > 100:
            #print('R greater than 100 on day',j)
        S_value = S(R_value,args)
        I_value = I(R_value,S_value,args)
        R_values.append(R_value)
        S_values.append(S_value)
        I_values.append(I_value)
        if (I_value > I_value_max):
            I_value_max = I_value
            I_value_max_day = j
    print('--- I value max =', I_value_max, 'on day', I_value_max_day)
    muy = args['beta']/args['gamma']
    I_value_max_true = ((-1) / (muy)) * math.log(muy) - (1/muy) + (args['N'] - args['n']) + args['n'] - (1/muy) * math.log(args['n'])
    plt.figure(figsize=(20,12))
    plt.plot(days,R_values,label='R', linewidth=4.5)
    plt.plot(days,S_values,label='S', linewidth=4.5)
    plt.plot(days,I_values,label='I', linewidth=4.5)
    plt.title('SIR model ('\
            + 'gamma: ' + str(args['gamma'])\
            + ', beta: ' + str(args['beta'])\
            + ')', fontsize=25)
    plt.xlabel('Days', fontsize=25)
    plt.ylabel('People', fontsize=25)
    plt.tick_params(labelsize=25)
    plt.legend(prop = { "size": 40 }, loc='center right')
    plt.show()
    plt.savefig('plot_method3.png')
    print('--- plotted ' + 'P: ' + str(args['P']) + ', K: ' + str(args['K']) + ' ---') 
    print('--- I value max =', I_value_max, 'on day', I_value_max_day)
    print('--- I value max true =', I_value_max_true)
    plt.clf()


if __name__=='__main__':
    #start_time = time.time()
    main()
    #stop_time = time.time()
    #execution_time = stop_time - start_time
    #print('Execution time: ',execution_time)
