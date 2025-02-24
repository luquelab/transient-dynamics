#This script generates inputs that are functions of the basic parameters
import math

r=0.9
d=3e-8
c=150
m=0.0028
tau=140
K=0

weight_r=r*tau
weight_m=m*tau
tau_r=1/r
tau_m=1/m
t_f=tau

B_c=1/(c*d*tau)
P_c=1/(d*tau)

Beq=m/(c*d);Peq=(r/d)
try:
        Peq_log=(r/d)*(1-m/(K*c*d))
except ZeroDivisionError:
        print("no carrying capacity")


Outputs=[r,d,c,m,tau,weight_r,weight_m,tau_r,tau_m,B_c,P_c,Beq,Peq]

with open('Derived_Parameters.txt','w', encoding="utf-8") as parameters_to_save:   
        for item in Outputs:
            parameters_to_save.write( str(item) + "\n")


