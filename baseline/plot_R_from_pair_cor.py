import matplotlib.pyplot as plt
import numpy as np

#jedine dvije varijable koje se mjenjaju po potrebi
Tc = 0.2
dim = 2

def find_zero_crossings(filename):
    data = np.loadtxt(filename, delimiter=' ')
    x = data[:, 0]
    y = data[:, 1]

    # Pronađi indekse gdje se događa promjena predznaka
    sign_changes = np.where(np.diff(np.sign(y)) != 0)[0]

    # Izračunaj nultočke pomoću linearnog pristupa između točaka
    zero_crossings = []
    for i in sign_changes:
        x0, x1 = x[i], x[i + 1]
        y0, y1 = y[i], y[i + 1]
        
        # Interpolacija nultih točaka
        x_zero = x0 - y0 * (x1 - x0) / (y1 - y0)
        zero_crossings.append(x_zero)

    return zero_crossings[0]
    
    # Plotaj graf
    """
    plt.plot(x, y, label='Funkcija')
    plt.scatter(zero_crossings, np.zeros_like(zero_crossings), color='red', label='Nultočke')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Graf funkcije s označenim nultočkama')
    plt.show()
    """
    #print("Nultočke su na x = ", zero_crossings)

R = []
t = []
t_exp = []

power_of_two = [2**n for n in range(34)]
thousand_times_n = [9000+1000*n for n in range(34)]
list = []
for i in range(34):
    list.append(power_of_two[i])
    list.append(thousand_times_n[i])

for i in list:
    try:
        R.append(find_zero_crossings( f"Data_{dim}d-T={Tc}Tc/vector_{i}.txt"))
        t.append(pow(i,(1/3)))
        t_exp.append(i)
    except:
        print(f"Datoteka: Data_{dim}d-T={Tc}Tc/vector_{i}.txt nije pronadena")

#Metoda najmanjih kvadrata
coefficients = np.polyfit(t, R, 1)
slope, intercept = coefficients
def line(x):
    return slope * x + intercept
x_fit = np.linspace(min(t), max(t), 100)
y_fit = line(x_fit)

plt.scatter(t,R,color='red',label =' Veličina domene kao nultočka C(r)')
plt.plot(x_fit, y_fit, color='blue', label='Pravac - najmanji kvadrati')
plt.title("Promjena veličine domene sustava vs $t^{1/3}$ \n izračunate kao nultočka parne korelacijske funkcije")
plt.ylabel("R - velicina domene")
plt.xlabel("$t[num_{mcs}]^{1/3}$")
plt.legend(loc = 'upper left')
plt.savefig( f"Data_{dim}d-T={Tc}Tc/R_vs_t/R-vs-$t^3$(correlation)_{dim}d.png")
plt.close()

coefficients = np.polyfit(np.log(t_exp), np.log(R), 1)
slope, intercept = coefficients
def line(x):
    return slope * x + intercept
x_fit = np.linspace(0, max(np.log(t_exp)), 100)
y_fit = line(x_fit)

coefficients2 = np.polyfit(np.log(t_exp[21:]), np.log(R[21:]), 1)
slope2, intercept2 = coefficients2
def line(x):
    return slope2 * x + intercept2
x_fit2 = np.linspace(min(np.log(t_exp[21:])), max(np.log(t_exp)), 100)
y_fit2 = line(x_fit2)

plt.plot(x_fit, y_fit, color='blue', label=fr'Najmanji kvadrati: y = {slope:.2f}$\cdot$ x {intercept:+.2f}')
plt.plot(x_fit2, y_fit2, color='red', label=fr'Najmanji kvadrati: y = {slope2:.2f}$\cdot$ x {intercept2:+.2f}')
plt.scatter(np.log(t_exp), np.log(R), marker='o', linestyle='-', label = 'Veličina domene')
plt.xlabel('log(t[$num_{mcs}$])')
plt.ylabel('log(R)')
plt.legend(loc = "upper left")
plt.title("log(t)-log(R) plot promjene veličine domene sustava \n izračunate kao nultočka parne korelacijske funkcije")
plt.savefig( f"Data_{dim}d-T={Tc}Tc/R_vs_t/log(R)-vs-log(t)(correlation)_{dim}d.png")

print(f"Eksponenet log-log je {slope}")