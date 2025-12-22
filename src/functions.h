using namespace std;

void printLattice(vector<vector<int>>& lattice) {
    int L = lattice.size();  // Mreža kvadratna, tj. LxL
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            cout << lattice[i][j] << " ";  
        }
        cout << endl;
    }
}
vector<vector<int>> createOrderedLattice(int num_A, int num_B, int L) {
    // Provjeri da li brojevi odgovaraju veličini matrice
    if (num_A + num_B > L * L) {
        throw invalid_argument("Broj elemenata num_A i num_B ne može biti veći od veličine matrice.");
    }

    // Inicijaliziraj praznu matricu (lattice) sa 0
    vector<vector<int>> lattice(L, vector<int>(L, 0));

    // Redom popuni matricu sa num_A i num_B
    int index = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (index < num_A) {
                lattice[i][j] = 1;  // Postavi 'A' (1)
            } else if (index < num_A + num_B) {
                lattice[i][j] = -1; // Postavi 'B' (-1)
            }
            index++;
        }
    }

    return lattice;
}


vector<vector<int>> createRandomLattice(int num_A, int num_B, int L) {
    // Stvori niz sa 'A' i 'B' vrijednostima
    vector<int> array;
    array.insert(array.end(), num_A, 1);  // 'A' -> 1
    array.insert(array.end(), num_B, -1); // 'B' -> -1

    // Postavi seed za generator slučajnih brojeva
    auto now = chrono::high_resolution_clock::now();
    auto seed = chrono::duration_cast<chrono::milliseconds>(now.time_since_epoch()).count();
    
    // Kreirajte generator slučajnih brojeva i distribuciju
    mt19937 g(seed);
    uniform_int_distribution<int> dist(0, L-1);  // Distribucija za brojeve između 0 i 60

    // Pomiješaj elemente niza
    shuffle(array.begin(), array.end(), g);

    // Pretvori u 2D niz (lattice) 
    vector<vector<int>> lattice(L, vector<int>(L));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            lattice[i][j] = array[i * L + j];
        }
    }
    // Dodaj 'V' (0) na nasumično mjesto
    
    int random_i = dist(g);
    int random_j = dist(g); 
    lattice[random_i][random_j] = 0;

    return lattice;
};

void saveLatticeToFile(const vector<vector<int>>& lattice, const string& filename) {
    ofstream outFile(filename);
    
    if (outFile.is_open()) {
        for (const auto& row : lattice) {
            for (int elem : row) {
                outFile << elem << " ";
            }
            outFile << "\n";  // Novi redak nakon svakog retka matrice
        }
        outFile.close();
        cout << "Lattice saved to " << filename << endl;
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }
};

int getEnergy(const vector<vector<int>>& lattice, int L, int J) {
    int energy = 0;
    vector<pair<int, int>> neighbor = {{1, 0}, {0, 1}};  // Mogući potezi (desno, dolje)

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (lattice[i][j] != 0) {  // Ako nije 'V' (0)
                for (auto move : neighbor) {
                    int ni = (i + move.first);
                    int nj = (j + move.second);

                    //Primijeni periodične rubne uvjete
                    
                    ni = ni % L;
                    nj = nj % L;

                    if (lattice[i][j] == -1) {  // 'A' ekvivalent
                        if (lattice[ni][nj] == -1) {
                            energy -= 1;
                        } else if (lattice[ni][nj] == 1) {  // 'B' ekvivalent
                            energy += 1;
                        }
                    } 
                    else if (lattice[i][j] == 1) {  // 'B' ekvivalent
                        if (lattice[ni][nj] == -1) {
                            energy += 1;
                        } else if (lattice[ni][nj] == 1) {
                            energy -= 1;
                        }
                    }
                }
            }
        }
    }

    return energy * J;
}

void apply_periodic_boundary_conditions(pair<int, int>& pos, int L) {
    int x = (pos.first + L) % L; 
    int y = (pos.second + L) % L;
    pos = {x,y};
}

void swapLatticePositions(vector<vector<int>>& lattice, pair<int, int>& vacancy_pos,pair<int, int>& new_pos) {
    int temp = lattice[vacancy_pos.first][vacancy_pos.second];
    
    // Swap vrijednosti
    lattice[vacancy_pos.first][vacancy_pos.second] = lattice[new_pos.first][new_pos.second];
    lattice[new_pos.first][new_pos.second] = temp;
}

void find_vacancy(vector<vector<int>>& lattice, int& L,pair<int, int>& vacancy){
    
    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            if (lattice[x][y] == 0) {
                vacancy.first = x;
                vacancy.second = y;};
            }
        }
    }

vector<double> pair_correlation_function(const vector<vector<int>>& lattice, int L) {
    int max_ind = L / 2;  // Maksimalni indeks za udaljenosti je L/2
    vector<double> f(max_ind + 1, 0.0);  // Inicijaliziraj polje korelacijske funkcije s 0.0

    // Glavna petlja koja prolazi kroz sve parove u rešetki
    for (int ir = 1; ir < L; ++ir) {
        for (int jr = 1; jr < L; ++jr) {
            for (int k = 0; k < L; ++k) {
                if (k != ir) {
                    int ind;
                    if ((k - ir) >= L / 2) {
                        ind = abs(k - ir - L);
                    } else if ((k - ir) <= -L / 2) {
                        ind = k - ir + L;
                    } else {
                        ind = abs(k - ir);
                    }
                    f[ind] += lattice[ir][jr] * lattice[k][jr];
                }

                if (k != jr) {
                    int ind;
                    if ((k - jr) >= L / 2) {
                        ind = abs(k - jr - L);
                    } else if ((k - jr) <= -L / 2) {
                        ind = k - jr + L;
                    } else {
                        ind = abs(k - jr);
                    }
                    f[ind] += lattice[ir][jr] * lattice[ir][k];
                }
            }
        }
    }

    // Normiranje vektora f
    double sum_f = accumulate(f.begin(), f.end(), 0.0);  // Izračunaj sumu svih elemenata u f
    if (sum_f != 0.0) {  // Izbjegni dijeljenje s nulom
        for (auto& val : f) {
            val /= sum_f;  // Podijeli svaki element s ukupnom sumom
        }
    }

    return f;
}

void calculate_R(vector<double>& energies,int& L, vector<double>& results){
    for (int value : energies) {
        results.push_back(2/((value / pow(L,2))+2)); 
    } 
}


void writeVectorsToFile(const vector<int>& vec1, const vector<double>& vec2, string filename) {
    // Provjeri da li su vektori iste veličine
    if (vec1.size() != vec2.size()) {
        cerr << "Vektori moraju biti iste veličine!" << endl;
        return; // Izlazi iz funkcije s greškom
    }

    // Otvori datoteku za pisanje
    ofstream outfile(filename);

    // Provjeri je li datoteka uspješno otvorena
    if (!outfile.is_open()) {
        cerr << "Ne mogu otvoriti datoteku za pisanje!" << endl;
        return; // Izlazi iz funkcije s greškom
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        outfile << vec1[i] << " " << vec2[i] << endl;
    }

    outfile.close();

    cout << "Podaci su uspjesno zapisani u " << filename << endl;
}

void saveCorrelationsToFiles(const vector<vector<double>>& pair_correlation,vector<int>& time,float factor,int br_decimala) {
    for (size_t i = 0; i < pair_correlation.size(); ++i) {
        char filename[100];
        sprintf(filename, "Data_2d-T=%.*fTc/vector_%d.txt", br_decimala, factor,time[i]);  
        //string filename = "Data-T=" + to_string(factor) + "Tc/" + "vector_" + to_string(time[i]) + ".txt";
        ofstream file(filename);

        if (!file) {
            cerr << "Ne mogu otvoriti datoteku " << filename << endl;
            continue;
        }

        const vector<double>& vec = pair_correlation[i];

        // Zapisivanje indeksa i vrijednosti u datoteku
        for (size_t j = 1; j < vec.size(); ++j) {
            file << j << " " << fixed << setprecision(2) << vec[j] << endl;
        }

        file.close();
        cout << "Zapisana korelacija za t:"<< time[i] <<endl;
    }
}