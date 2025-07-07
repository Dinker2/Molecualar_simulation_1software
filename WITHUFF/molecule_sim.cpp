#include <bits/stdc++.h>
using namespace std;

struct Atom {
    string symbol;
    double x, y, z;
    double electronegativity;
};

map<string, double> covalent_radius = {
    {"H", 0.31}, {"C", 0.76}, {"N", 0.71}, {"O", 0.66},
    {"F", 0.57}, {"P", 1.07}, {"S", 1.05}, {"Cl", 1.02}
};

map<string, double> vdw_radius = {
    {"H", 1.20}, {"C", 1.70}, {"N", 1.55}, {"O", 1.52},
    {"F", 1.47}, {"P", 1.80}, {"S", 1.80}, {"Cl", 1.75}
};

map<string, int> valency = {
    {"H", 1}, {"C", 4}, {"N", 3}, {"O", 2},
    {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}
};

map<string, double> electronegativity = {
    {"H", 2.20}, {"C", 2.55}, {"N", 3.04}, {"O", 3.44},
    {"F", 3.98}, {"P", 2.19}, {"S", 2.58}, {"Cl", 3.16}
};

double distance(const Atom& a, const Atom& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

vector<Atom> translate(const vector<Atom>& unit, double dx, double dy, double dz) {
    vector<Atom> shifted;
    for (const auto& a : unit) {
        shifted.push_back({a.symbol, a.x + dx, a.y + dy, a.z + dz, a.electronegativity});
    }
    return shifted;
}

vector<Atom> random_rotate(const vector<Atom>& unit) {
    double ax = ((double)rand() / RAND_MAX) * 2 * M_PI;
    double ay = ((double)rand() / RAND_MAX) * 2 * M_PI;
    double az = ((double)rand() / RAND_MAX) * 2 * M_PI;

    double cx = cos(ax), sx = sin(ax);
    double cy = cos(ay), sy = sin(ay);
    double cz = cos(az), sz = sin(az);

    vector<Atom> rotated;
    for (auto a : unit) {
        double x = a.x, y = a.y, z = a.z;

        double y1 = cy * y - sy * z;
        double z1 = sy * y + cy * z;

        double x2 = cx * x + sx * z1;
        double z2 = -sx * x + cx * z1;

        double x3 = cz * x2 - sz * y1;
        double y3 = sz * x2 + cz * y1;

        rotated.push_back({a.symbol, x3, y3, z2, a.electronegativity});
    }

    return rotated;
}

bool is_valid(const vector<Atom>& new_unit, const vector<Atom>& existing, double min_dist) {
    for (auto& a : new_unit) {
        for (auto& b : existing) {
            if (distance(a, b) < min_dist) return false;
        }
    }
    return true;
}

vector<pair<int, int>> infer_bonds(const vector<Atom>& atoms, double tol = 0.45) {
    vector<pair<int, int>> bonds;
    for (int i = 0; i < atoms.size(); ++i) {
        for (int j = i + 1; j < atoms.size(); ++j) {
            double cutoff = covalent_radius[atoms[i].symbol] + covalent_radius[atoms[j].symbol] + tol;
            if (distance(atoms[i], atoms[j]) <= cutoff) {
                bonds.emplace_back(i, j);
            }
        }
    }
    return bonds;
}

int main() {
    srand(time(0));

    // === Benzene Template ===
    vector<Atom> reference = {
   {"O" , 0.00 , -0.064 , 0.00},
  {"H" , 0.816 , 0.513 , 0.00},
  {"H" , -0.816 , 0.513 , 0.00}
    };

    for (auto& a : reference) {
        a.electronegativity = electronegativity.count(a.symbol) ? electronegativity[a.symbol] : 2.5;
    }

    // === Simulation Settings ===
    int no_of_units = 33;
    double cube_size = 10.0;         // Ångström
    double min_dist = 2.5;            // Ångström

    // === New Input: Density in g/cm³ and Molecular Weight ===
    double density, mol_weight;
    /*Enter material density (in g/cm^3)*/
    density= 1;
    /*Enter molecular weight of one molecule (in g/mol)*/
    mol_weight= 18;

    // === Compute max units that can fit ===
    const double AVOGADRO = 6.022e23;
    double box_volume_cm3 = (cube_size * cube_size * cube_size) * 1e-24; // Å³ to cm³
    double total_mass = density * box_volume_cm3; // grams
    int max_possible_units = static_cast<int>((total_mass * AVOGADRO) / mol_weight);

    cout << "Estimated maximum number of molecules that can fit in the box: " << max_possible_units << "\n";

    // === Molecule Placement ===
    vector<Atom> all_atoms;
    int tries = 0, max_attempts = 10000;
    int placed_units = 0;

    while (placed_units < no_of_units && tries < max_attempts) {
        double dx = ((double)rand() / RAND_MAX) * cube_size;
        double dy = ((double)rand() / RAND_MAX) * cube_size;
        double dz = ((double)rand() / RAND_MAX) * cube_size;

        vector<Atom> unit = translate(random_rotate(reference), dx, dy, dz);

        if (is_valid(unit, all_atoms, min_dist)) {
            all_atoms.insert(all_atoms.end(), unit.begin(), unit.end());
            placed_units++;
        }

        ++tries;
    }

    // === Output to file ===
    ofstream fout("points_with_symbols.txt");
    fout << fixed << setprecision(3);
    for (auto& a : all_atoms) {
        fout << a.symbol << " "
             << a.x << " " << a.y << " " << a.z << " "
             << a.electronegativity << "\n";
    }
    fout.close();

    cout << "Simulation complete.\n";
    cout << "Placed units: " << placed_units << "\n";
    cout << "Output (coordinates in Å, EN in Pauling) written to points_with_symbols.txt\n";

    return 0;
}
