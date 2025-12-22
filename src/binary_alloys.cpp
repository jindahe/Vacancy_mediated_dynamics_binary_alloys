#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cmath>
#include <cstdlib>
#include "functions.h"
#include <algorithm>  // Za std::shuffle
#include <random>     // Za generator slučajnih brojeva
#include <ctime>      // Za seeding random generatora
#include <vector>
#include <chrono>
#include <random>
#include <cstdio>

using namespace std;

void move(vector<vector<int>>& lattice, int& start_energy, double J, double T, int L,int num_mc, int num_steps,vector<double>& energies_in_time,vector<int>& times,vector<vector<double>>& pair_correlation,vector<int>& time_correlation,vector<int>& powers_of_two,float factor,int br_decimala) {
    double random_number,dE;
    int x,y,new_energy,old_energy,sum_energy,random_index;
    vector<pair<int, int>> neighbor = {{-1,0},{1, 0},{0,-1}, {0, 1}};  
    pair<int, int> random_move, vacancy_pos, new_pos;
    vacancy_pos = {0,0};

    sum_energy = start_energy;
    auto now = chrono::high_resolution_clock::now();
    auto seed = chrono::duration_cast<chrono::nanoseconds>(now.time_since_epoch()).count();
    random_device rd;
    seed ^= rd();  // XOR sa nasumičnim brojem

    mt19937 g(seed);
    uniform_int_distribution<int> dist(0, neighbor.size() - 1);
    // Inicijaliziraj generator slučajnih brojeva sa seed-om trenutnog vremena
    mt19937 gen(seed);
    uniform_real_distribution<double> dis(0.0, 1.0); // Uniformna distribucija između 0 i 1
    
    find_vacancy(lattice,L,vacancy_pos); 

    for(int i = 0;i < num_mc; i++){
        for(int j = 0;j < num_steps; j++){   
            
            random_index = dist(g);
            random_move = neighbor[random_index];
                 
            x = random_move.first;
            y = random_move.second;

            new_pos = {0,0}; 
            new_pos.first = vacancy_pos.first + x; 
            new_pos.second = vacancy_pos.second + y;
            if(new_pos.first < 0 || new_pos.first > (L-1) ||new_pos.second < 0 || new_pos.second > (L-1)){
                apply_periodic_boundary_conditions(new_pos,L);
            }

            old_energy = -J * (
            lattice[new_pos.first][new_pos.second] * lattice[(new_pos.first - 1 + L) % L][new_pos.second] +
            lattice[new_pos.first][new_pos.second] * lattice[(new_pos.first + 1) % L][new_pos.second] +
            lattice[new_pos.first][new_pos.second] * lattice[new_pos.first][(new_pos.second + 1) % L] +
            lattice[new_pos.first][new_pos.second] * lattice[new_pos.first][(new_pos.second - 1 + L) % L]
            );

            swapLatticePositions(lattice, vacancy_pos, new_pos);
            new_energy = -J * (
            lattice[vacancy_pos.first][vacancy_pos.second] * lattice[(vacancy_pos.first - 1 + L) % L][vacancy_pos.second] +
            lattice[vacancy_pos.first][vacancy_pos.second] * lattice[(vacancy_pos.first + 1) % L][vacancy_pos.second] +
            lattice[vacancy_pos.first][vacancy_pos.second] * lattice[vacancy_pos.first][(vacancy_pos.second + 1) % L] +
            lattice[vacancy_pos.first][vacancy_pos.second] * lattice[vacancy_pos.first][(vacancy_pos.second - 1 + L) % L]
            );  

            dE = new_energy - old_energy;
            
            if(dE<=0){
                start_energy += dE;
                vacancy_pos.first = new_pos.first;
                vacancy_pos.second = new_pos.second;
            }
            else if (dis(gen) <= exp(-dE/T))
            {
                start_energy += dE;
                vacancy_pos.first = new_pos.first;
                vacancy_pos.second = new_pos.second;
            }
            else {swapLatticePositions(lattice, vacancy_pos, new_pos);
            };
        }
        if (find(powers_of_two.begin(), powers_of_two.end(), i) != powers_of_two.end()) {
            char filename[100];
            sprintf(filename, "Data_2d-T=%.*fTc/lattice_output_t-%d.txt", br_decimala, factor,i);
            saveLatticeToFile(lattice, filename);
        }
        if(i %1000 ==0){
            cout<<"Broj msc:"<< i <<endl;
            //cout << getEnergy(lattice, L, J) - start_energy << endl; //provjera izracuna energije
        }
        if (find(powers_of_two.begin(), powers_of_two.end(), i) != powers_of_two.end() || (i > 8000 && i % 1000 == 0)) {
            pair_correlation.push_back(pair_correlation_function(lattice, L));
            time_correlation.push_back(i);
        }
        sum_energy += start_energy;
        energies_in_time.push_back(sum_energy/(i+1));
        times.push_back(i+1);
    }
}


int main(){
    int L = 128; //velicina resetke
    int N = L*L;
    float Tc = 2.26; //kriticna temperatura
    float factor = 0.5;
    int br_decimala = 1; //broj decimala u factor kako bi se ispravno otvorila mapa za zapisivanje
    double T = Tc*factor; //trenutna temperatura
    int num_A = N / 2, num_B = N / 2; //broj atoma A i B
    int num_mc = 33000; //broj mc koraka
    int num_steps = L*L; //broj koraka u jednom mc
    int J = 1; //koeficijent energije interakcije 
    
    vector<int> powers_of_two; //trenutci u kojima zapisujemo konfiguraciju
    for (int n = 0; n < 30; ++n) {  
            int value = pow(2, n);
            powers_of_two.push_back(value);
        }

    //Vectori za zapisivanje rezultata simmulacije
    vector<double> R,energies_in_time; 
    //vector<vector<int>> lattice = createOrderedLattice(num_A, num_B, L); //testiranje ispravnosti racuna energije
    vector<vector<double>> pair_correlation;
    vector<int> time_correlation,times;

    vector<vector<int>> lattice = createRandomLattice(num_A, num_B, L);

    char filename[100];
    sprintf(filename, "Data_2d-T=%.*fTc/lattice_output_t-0.txt", br_decimala, factor);  
    saveLatticeToFile(lattice, filename);

    int start_energy = getEnergy(lattice, L, J);
    cout << "Total energy: " << start_energy << endl;
    
    move(lattice, start_energy, J, T, L, num_mc,num_steps,energies_in_time,times,pair_correlation,time_correlation,powers_of_two,factor,br_decimala);
    calculate_R(energies_in_time,L,R);
    
    sprintf(filename, "Data_2d-T=%.*fTc/t_vs_R.txt",br_decimala,factor); 
    writeVectorsToFile(times, R, filename);
    saveCorrelationsToFiles(pair_correlation,time_correlation,factor,br_decimala);
}
