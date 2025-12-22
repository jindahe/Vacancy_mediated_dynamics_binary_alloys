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
vector<vector<vector<int>>> createOrderedLattice(int num_A, int num_B, int L) {
    // Provjeri da li brojevi odgovaraju veličini matrice
    if (num_A + num_B > L * L * L) {
        throw invalid_argument("Broj elemenata num_A i num_B ne može biti veći od veličine matrice.");
    }

    // Inicijaliziraj praznu matricu (lattice) sa 0
    vector<vector<vector<int>>> lattice(L, vector<vector<int>>(L, vector<int>(L, 0)));

    // Redom popuni matricu sa num_A i num_B
    int index = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                if (index < num_A) {
                    lattice[i][j][k] = 1;  // Postavi 'A' (1)
                } else if (index < num_A + num_B) {
                    lattice[i][j][k] = -1; // Postavi 'B' (-1)
                }
                index++;
            }
        }
    }

    return lattice;
}


vector<vector<vector<int>>> createRandomLattice(int num_A, int num_B, int L) {
    // Stvori niz sa 'A' i 'B' vrijednostima
    vector<int> array;
    array.insert(array.end(), num_A, 1);  // 'A' -> 1
    array.insert(array.end(), num_B, -1); // 'B' -> -1

    // Postavi seed za generator slučajnih brojeva
    auto now = chrono::high_resolution_clock::now();
    auto seed = chrono::duration_cast<chrono::milliseconds>(now.time_since_epoch()).count();
    
    // Kreirajte generator slučajnih brojeva i distribuciju
    mt19937 g(seed);
    uniform_int_distribution<int> dist(0, L-1);  // Distribucija za brojeve između 0 i L-1

    // Pomiješaj elemente niza
    shuffle(array.begin(), array.end(), g);

    // Pretvori u 3D niz (lattice) 
    vector<vector<vector<int>>> lattice(L, vector<vector<int>>(L, vector<int>(L)));
    int index = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                lattice[i][j][k] = array[index++];
            }
        }
    }

    // Dodaj 'V' (0) na nasumično mjesto
    int random_i = dist(g);
    int random_j = dist(g);
    int random_k = dist(g);
    lattice[random_i][random_j][random_k] = 0;

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

int getEnergy(const vector<vector<vector<int>>>& lattice, int L, int J) {
    int energy = 0;
    vector<tuple<int, int, int>> neighbor = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};  // Mogući potezi (desno, dolje, naprijed)

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                if (lattice[i][j][k] != 0) {  // Ako nije 'V' (0)
                    for (auto move : neighbor) {
                        int ni = (i + get<0>(move)) % L;
                        int nj = (j + get<1>(move)) % L;
                        int nk = (k + get<2>(move)) % L;

                        if (lattice[i][j][k] == -1) {  // 'A' ekvivalent
                            if (lattice[ni][nj][nk] == -1) {
                                energy -= 1;
                            } else if (lattice[ni][nj][nk] == 1) {  // 'B' ekvivalent
                                energy += 1;
                            }
                        } else if (lattice[i][j][k] == 1) {  // 'B' ekvivalent
                            if (lattice[ni][nj][nk] == -1) {
                                energy += 1;
                            } else if (lattice[ni][nj][nk] == 1) {
                                energy -= 1;
                            }
                        }
                    }
                }
            }
        }
    }

    return energy * J;
}


void apply_periodic_boundary_conditions(tuple<int, int, int>& pos, int L) {
    int x = (get<0>(pos) + L) % L;
    int y = (get<1>(pos) + L) % L;
    int z = (get<2>(pos) + L) % L;
    pos = {x, y, z};
}


void swapLatticePositions(vector<vector<vector<int>>>& lattice, tuple<int, int, int>& vacancy_pos, tuple<int, int, int>& new_pos) {
    int temp = lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)];
    
    // Swap vrijednosti
    lattice[get<0>(vacancy_pos)][get<1>(vacancy_pos)][get<2>(vacancy_pos)] = lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)];
    lattice[get<0>(new_pos)][get<1>(new_pos)][get<2>(new_pos)] = temp;
}


void find_vacancy(vector<vector<vector<int>>>& lattice, int& L, tuple<int, int, int>& vacancy) {
    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            for (int z = 0; z < L; ++z) {
                if (lattice[x][y][z] == 0) {
                    get<0>(vacancy) = x;
                    get<1>(vacancy) = y;
                    get<2>(vacancy) = z;
                    return;  // Exit once the vacancy is found
                }
            }
        }
    }
}


vector<double> pair_correlation_function(const vector<vector<vector<int>>>& lattice, int L) {
    int max_ind = L / 2;  // Maksimalni indeks za udaljenosti je L/2
    vector<double> f(max_ind + 1, 0.0);  // Inicijaliziraj polje korelacijske funkcije s 0.0

    // Glavna petlja koja prolazi kroz sve parove u rešetki
    for (int ir = 1; ir < L; ++ir) {
        for (int jr = 1; jr < L; ++jr) {
            for(int zr=1; zr < L;++zr){
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
                    f[ind] += lattice[ir][jr][zr] * lattice[k][jr][zr];
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
                    f[ind] += lattice[ir][jr][zr] * lattice[ir][k][zr];
                }
                if(k != zr){
                    int ind;
                    if ((k - zr) >= L / 2) {
                        ind = abs(k - zr - L);
                    } else if ((k - zr) <= -L / 2) {
                        ind = k - zr + L;
                    } else {
                        ind = abs(k - zr);
                    }
                    f[ind] += lattice[ir][jr][zr] * lattice[ir][jr][k];
                }
            }
            }
        }
    }

    return f;
}


void calculate_R(vector<double>& energies,int& L, vector<double>& results){
    for (int value : energies) {
        results.push_back(3/((value / pow(L,3))+3)); 
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

    //cout << "Podaci su uspjesno zapisani u " << filename << endl;
}

void saveCorrelationsToFiles(const vector<vector<double>>& pair_correlation,vector<int>& time,int br_decimala, float factor) {
    for (size_t i = 0; i < pair_correlation.size(); ++i) {
        char filename[100];
        sprintf(filename, "Data_3d-T=%.*fTc/vector_%d.txt", br_decimala, factor,time[i]);
        //string filename = "vector_" + to_string(time[i]) + ".txt";
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
        //cout << "Zapisana korelacija za t:"<< time[i] <<endl;
    }
}