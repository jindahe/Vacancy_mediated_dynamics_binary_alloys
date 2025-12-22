import matplotlib.pyplot as plt

#jedine dvije varijable koje se mjenjaju po potrebi
Tc = 0.5
dim = 2

power_of_two = [2**n for n in range(34)]
thousand_times_n = [9000+1000*n for n in range(34)]
list = []
for i in range(34):
    list.append(power_of_two[i])
    list.append(thousand_times_n[i])
    
for i in list:
    filename = f"Data_{dim}d-T={Tc}Tc/vector_{i}.txt"  

    vec1 = []
    vec2 = []
    j = 0
    try:
        # Otvorite datoteku za čitanje
        with open(filename, 'r') as file:
            for line in file:
                # Razdvoji redak na dva dijela po razmaku
                if j>0:
                    parts = line.split()
                    if len(parts) == 2:
                        # Dodaj podatke u odgovarajuće vektore
                        vec1.append(float(parts[0]))
                        vec2.append(float(parts[1]))
                    else:
                        print(f"Upozorenje: Neregularan redak u datoteci: {line.strip()}")
                j += 1
    except FileNotFoundError:
        print(f"Datoteka {filename} nije pronađena!")
        continue
    except Exception as e:
        print(f"Dogodila se greška: {e}")
        continue

    plt.plot(vec1,vec2)
    plt.title(f"Parna korelacijska funkcija za t={i} mcs")
    plt.ylabel("C(r)")
    plt.xlabel("r[Atomska udaljenost u rešetci]")
    plt.savefig(f"Data_{dim}d-T={Tc}Tc/C(r)/Pair_correlation-{i}.png")
    plt.close()
