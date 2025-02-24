import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def f(xs, t, ps):
    """Lotka-Volterra predator-prey model."""

    try:
        a = ps['a'].value
        b = ps['b'].value
        c = ps['c'].value
        d = ps['d'].value
    except:
        a, b, c, d = ps

    x, y = xs
    return [a*x - b*x*y, c*x*y - d*y]

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """

    x = odeint(f, x0, t, args=(ps,))


    return x


def residual(ps, ts, data):
    x0 = ps['x0'].value, ps['y0'].value
    model = g(ts, x0, ps)

    return (model - data).ravel()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Time and initial conditions
#============================
Length=100
t = np.linspace(0, 10, Length)
x0 = np.array([1,1])
#============================

#Parameter values
#==========================================
a, b, c, d = 3,1,1,1
true_params = np.array((a, b, c, d))

#Solve model
model = g(t, x0, true_params)
#Add noise
noise = np.random.normal(size=model.shape)
#data=model+noise
data=model

#Small data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Time_indices=np.linspace(0,Length-1,Length)
Random_Indices=np.random.choice(Time_indices, 100, replace=False)
Random_Indices=[int(random_index) for random_index in Random_Indices]
Random_Indices=np.sort(Random_Indices)
Small_data=data[Random_Indices]
Small_time=t[Random_Indices]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#==========================================

# set parameters incluing bounds
#===========================================================
params = Parameters()
params.add('x0', value= float(data[0, 0]), min=0, max=10)
params.add('y0', value=float(data[0, 1]), min=0, max=10)
params.add('a', value=2.0, min=0, max=10)
params.add('b', value=1.0, min=0, max=10)
params.add('c', value=1.0, min=0, max=10)
params.add('d', value=1.0, min=0, max=10)
#===========================================================


# fit model and find predicted values
#=====================================================================
# result = minimize(residual, params, args=(t, data), method='leastsq')
# final = data + result.residual.reshape(data.shape)

result=minimize(residual,params,args=(Small_time, Small_data),method='leastsq')
final = Small_data + result.residual.reshape(Small_data.shape)
#=====================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot data and fitted curves
# plt.plot(t, data, 'o')
# plt.plot(t, final, '-', linewidth=2);


plt.plot(t, model, '--')
plt.plot(Small_time, Small_data, 'o')
plt.plot(Small_time, final, '-', linewidth=2);

plt.show()
# display fitted statistics
report_fit(result)
print(result.nvarys)
print(result.ndata)
print(result.residual)
print(result.nfree)
print(result.chisqr)
print(result.redchi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
