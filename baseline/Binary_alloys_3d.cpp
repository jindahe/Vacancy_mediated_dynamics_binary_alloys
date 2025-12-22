#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <cmath>
#include <cstdlib>
#include "functions_3d.h"
#include <algorithm>  // Za std::shuffle
#include <random>     // Za generator slučajnih brojeva
#include <ctime>      // Za seeding random generatora
#include <vector>
#include <chrono>
#include <random>
#include <cstdio>
#include <tuple>

using namespace std;

void move(vector<vector<vector<int>>>& lattice, int& start_energy, double J, double T, int L,int num_mc, int num_steps,vector<double>& energies_in_time,vector<int>& times,vector<vector<double>>& pair_correlation,vector<int>& time_correlation,vector<int>& powers_of_two) {
    double random_number,dE;
    int x,y,z,new_energy,old_energy,random_index;
    long long int sum_energy;
    vector<tuple<int, int, int>> neighbor = {
        make_tuple(-1, 0, 0),
        make_tuple(1, 0, 0),
        make_tuple(0, -1, 0),
        make_tuple(0, 1, 0),
        make_tuple(0, 0, 1),
        make_tuple(0, 0, -1)
    };  
    tuple<int, int,int> random_move, vacancy_pos, new_pos;
    vacancy_pos = make_tuple(0,0,0);

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
    //cout << "tocka 1: "<<getEnergy(lattice, L, J) << endl;
    for(int i = 0;i < num_mc; i++){
        for(int j = 0;j < num_steps; j++){   
            
            random_index = dist(g);
            random_move = neighbor[random_index];
                 
            x = get<0>(random_move);
            y = get<1>(random_move);
            z = get<2>(random_move);

            new_pos =  make_tuple(0,0,0);
            get<0>(new_pos) = get<0>(vacancy_pos) + x;  
            get<1>(new_pos) = get<1>(vacancy_pos) + y; 
            get<2>(new_pos) = get<2>(vacancy_pos) + z;
           if (get<0>(new_pos) < 0 || get<0>(new_pos) > (L-1) || get<1>(new_pos) < 0 || get<1>(new_pos) > (L-1) || get<2>(new_pos) < 0 || get<2>(new_pos) > (L-1)) {
                apply_periodic_boundary_conditions(new_pos,L);
            }

            old_energy = -J * (
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[(get<0>(new_pos) - 1 + L) % L][get<1>(new_pos)][get<2>(new_pos)] +
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[(get<0>(new_pos) + 1) % L][get<1>(new_pos)][get<2>(new_pos)] +
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[get<0>(new_pos)][(get<1>(new_pos) + 1) % L][get<2>(new_pos)] +
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[get<0>(new_pos)][(get<1>(new_pos) - 1 + L) % L][get<2>(new_pos)] +
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[get<0>(new_pos)][get<1>(new_pos)][(get<2>(new_pos) + 1) % L] +
            lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] * lattice[get<0>(new_pos)][get<1>(new_pos)][(get<2>(new_pos) - 1 + L) % L]
        );
          
            swapLatticePositions(lattice, vacancy_pos, new_pos);

            new_energy = -J * (
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[(get<0>(vacancy_pos) - 1 + L) % L][get<1>(vacancy_pos)][get<2>(vacancy_pos)] +
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[(get<0>(vacancy_pos) + 1) % L][get<1>(vacancy_pos)][get<2>(vacancy_pos)] +
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[get<0>(vacancy_pos)][(get<1>(vacancy_pos) + 1) % L][get<2>(vacancy_pos)] +
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[get<0>(vacancy_pos)][(get<1>(vacancy_pos) - 1 + L) % L][get<2>(vacancy_pos)] +
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][(get<2>(vacancy_pos) + 1) % L] +
            lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] * lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][(get<2>(vacancy_pos) - 1 + L) % L]
            );


            if (dE <= 0) {
                start_energy += dE;
                get<0>(vacancy_pos) = get<0>(new_pos);
                get<1>(vacancy_pos) = get<1>(new_pos);
                get<2>(vacancy_pos) = get<2>(new_pos);
            } 
            else if (dis(gen) <= exp(-dE / T)) {
                start_energy += dE;
                get<0>(vacancy_pos) = get<0>(new_pos);
                get<1>(vacancy_pos) = get<1>(new_pos);
                get<2>(vacancy_pos) = get<2>(new_pos);
            }
            else {swapLatticePositions(lattice, vacancy_pos, new_pos);
            };
        }
        /*
        if (find(powers_of_two.begin(), powers_of_two.end(), i) != powers_of_two.end()) {
            char filename[100];
            sprintf(filename, "lattice_output_t-%d.txt", i);
            saveLatticeToFile(lattice, filename);
        }*/
        if(i %10 ==0){
            cout<<"Broj msc:"<< i <<endl;
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
    int L = 128;
    int N = L*L*L;
    float Tc = 2.26;
    float factor = 0.5;
    int br_decimala = 1;
    double T = Tc*factor;
    int num_A = N / 2, num_B = N / 2;
    int num_mc = 33000;
    int num_steps = L*L*L;
    int J = 1;
    vector<int> powers_of_two;
    for (int n = 0; n < 30; ++n) {  
            int value = pow(2, n);
            powers_of_two.push_back(value);
        }
    vector<double> R,energies_in_time;
    //vector<vector<int>> lattice = createOrderedLattice(num_A, num_B, L); //testiranje ispravnosti racuna energije
    vector<vector<vector<int>>> lattice = createRandomLattice(num_A, num_B, L);
    vector<vector<double>> pair_correlation;
    vector<int> time_correlation,times;
    
    //saveLatticeToFile(lattice, "lattice_output_t-0.txt");
    int start_energy = getEnergy(lattice, L, J);
    //cout << "Lattice size: " << L << endl;
    move(lattice, start_energy, J, T, L, num_mc,num_steps,energies_in_time,times,pair_correlation,time_correlation,powers_of_two);
    calculate_R(energies_in_time,L,R);
    
    char filename[100];
    sprintf(filename, "Data_3d-T=%.*fTc/t_vs_R.txt",br_decimala,factor); 
    writeVectorsToFile(times, R, filename);
    saveCorrelationsToFiles(pair_correlation,time_correlation,br_decimala,factor);
}
