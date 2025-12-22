import matplotlib.pyplot as plt
import numpy as np

#jedine dvije varijable koje se mjenjaju po potrebi
Tc = 0.5
dim = 2

filename = f"Data_{dim}d-T={Tc}Tc/t_vs_R.txt"

t_exp = []
t = []
R = []

try:
    with open(filename, 'r') as file:
        for line in file:
            # Razdvoji redak na dva dijela po razmaku
            parts = line.split()
            if len(parts) == 2:
                # Dodaj podatke u odgovarajuće vektore
                t_exp.append(float(parts[0]))
                t.append(float(parts[0])**(1/3))
                R.append(float(parts[1]))
            else:
                print(f"Upozorenje: Neregularan redak u datoteci: {line.strip()}")
except FileNotFoundError:
    print(f"Datoteka {filename} nije pronađena!")
except Exception as e:
    print(f"Dogodila se greška: {e}")

#Metoda najmanjih kvadrata
coefficients = np.polyfit(t, R, 1)
slope, intercept = coefficients
def line(x):
    return slope * x + intercept
x_fit = np.linspace(min(t), max(t), 100)
y_fit = line(x_fit)

plt.scatter(t,R,color='red',label =' Veličina domene')
plt.plot(x_fit, y_fit, color='blue', label='Pravac - najmanji kvadrati')
plt.title("Promjena veličine domene sustava vs $t^{1/3}$ \n izračunata preko ukupne energije sustava")
plt.ylabel("R - velicina domene")
plt.xlabel("$t[num_{mcs}]^{1/3}$")
plt.legend(loc = 'upper left')
plt.savefig(f"Data_{dim}d-T={Tc}Tc/R_vs_t/R-vs-$t^3$(energy)_{dim}d.png")
plt.close()

del t_exp[0]
del R[0]
coefficients = np.polyfit(np.log(t_exp), np.log(R), 1)
slope, intercept = coefficients
def line(x):
    return slope * x + intercept
x_fit = np.linspace(0, max(np.log(t_exp)), 100)
y_fit = line(x_fit)

plt.plot(x_fit, y_fit, color='blue', label=fr'Najmanji kvadrati: y = {slope:.2f}$\cdot$ x {intercept:+.2f}')
plt.scatter(np.log(t_exp), np.log(R), marker='o', linestyle='-', label = 'Veličina domene')
plt.legend(loc = 'upper left')
plt.title("log(t)-log(R) plot promjene veličine domene sustava \n izračunate preko ukupne energije sustava")
plt.xlabel('log(t[$num_{mcs}$])')
plt.ylabel('log(R)')
plt.savefig(f"Data_{dim}d-T={Tc}Tc/R_vs_t/log(R)-vs-log(t)(energy)_{dim}d.png")

print(f"Eksponenet log-log je {slope}")
