import matplotlib.pyplot as plt

#jedine dvije varijable koje se mjenjaju po potrebi
Tc = 0.5
dim = 2

list_i = [2**n for n in range (30)]
list_i.append(0)

for i in list_i: 
    lattice = []
    try:
        with open(f"Data_{dim}d-T={Tc}Tc/lattice_output_t-{i}.txt") as f:
            Lines = f.readlines()
            for line in Lines:
                lattice.append([])
                line = line.split(" ")
                for elem in line:
                    try:
                        br = int(elem)
                        lattice[-1].append(br)
                    except:
                        pass
                            
        plt.figure(figsize=(8, 8))
        plt.imshow(lattice, cmap='coolwarm', interpolation='none', vmin=-1, vmax=1)
        plt.title(f'State of the Lattice for t = {i}[mcs]\n(A atoms in blue, B atoms in red, Vacancy in white)')
        plt.colorbar(ticks=[-1, 0, 1], format=plt.FuncFormatter(lambda val, loc: ['V', "B","A"][int(val)]))
        plt.savefig(f"Data_{dim}d-T={Tc}Tc/Lattice_change_in_time/Izgled resetke za t ={i} mcs")
        print(f"Spremljen izgled resetke za t ={i} mcs")
        plt.close()
    except:
        pass
