nusp = 5436645 # Murilo Hiroaki Seko

##################################

import numpy as np
mnp = (nusp%1000)/1000

##################################
#   Questão 1
print("Questão 1")

c = round((1. + mnp)*1e-9, 13)
u = round((2. + mnp)*1e8, -5)
l = round((2. + mnp)*10, 2)

Z0 = 75. #50.
Rg = 27. #25.
Eg = 10.
Bl = 0.0000001

E = Eg*Z0/(Rg+Z0)
rhog = (Rg-Z0)/(Rg+Z0)

#       item a
print("item a)")

print("\t v1(t=0,2- us) =", round(E,5), "V")
print("\t v1(t=0,2+ us) =", round(-E*rhog,5), "V")
#print("\t v1(t=0,2+ us) =", round(E + E*(1 + rhog)*(1 - 2*np.exp((-0.2*1e-6 + 2*Bl)/(Z0*c))),5), "V")

#       item b
print()
print("item b)")

print("\t v1(t=inf) = v2(t=inf) =", round(Eg,5), "V")

#       item d
print()
print("item d)")

print("\t tau =", round(Z0*c*1e9,5), "ns")
print("\t v1(t=0,39999+ us) =", round(E + E*(1 + rhog)*(1 - 2*np.exp((-0.39999*1e-6 + 2*Bl)/(Z0*c))),5), "V")

#       item e
print()
print("item e)")

A = round((7. + mnp)*1e-3, 6) #round((5. + mnp)*1e-3, 6)

print("\t L =", round(Z0/u*1e9,5), "nH/m")
print("\t C =", round(1/(Z0*u)*1e12,5), "pF/m")
print("\t R =", round(A*Z0,5), "R/m")
print("\t G =", round(A/Z0*1e3,5), "mS/m")

# perdas = np.exp(-A*l)

print("\t v1(t=0,2001 us) =", round(E + E*(1 + rhog)*(1 - 2*np.exp((-0.2001*1e-6 + 2*Bl)/(Z0*c)))*np.exp(-A*2*l),5), "V")
print("\t v2(t=0,29999 us) =", round(2*E*(1 - np.exp((-0.29999*1e-6 + Bl)/(Z0*c)))*np.exp(-A*l),5), "V")


##################################
#   Questão 2
print()
print()
print("Questão 2")

Rl = round((4. + mnp)*100, 1)
u = round((2. + mnp)*1e8, -5)
l = round((2. + mnp)*10, 2)

Z0 = 75.
Rg = 5.
Eg = 10.
Bl = 0.0000001

# E = Eg*Z0/(Rg+Z0)
# rhog = (Rg-Z0)/(Rg+Z0)
# rhol = (Rl-Z0)/(Rl+Z0)

L = Z0/u
C = 1/Z0/u

print("\t C_10 =", round(C*l/10*1e12,5), "pF")
print("\t L_10 =", round(L*l/10*1e9,5), "nH")
print("\t C_20 =", round(C*l/20*1e12,5), "pF")
print("\t L_20 =", round(L*l/20*1e9,5), "nH")