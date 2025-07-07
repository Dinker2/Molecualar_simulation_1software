import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === Atom properties ===
atom_colors = {
    "H": "blue",
    "O": "red",
    "N": "white",
    "C": "red",
    "F": "green",
    "Cl": "green",
    "S": "yellow",
    "P": "orange"
}

atom_sizes = {
    "H": 100,
    "C": 100,
    "N": 250,
    "O": 100,
    "F": 240,
    "Cl": 240,
    "S": 280,
    "P": 280
}

# === Load atom data ===
atoms = []
with open("points_with_symbols.txt") as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbol = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((symbol, x, y, z))

# === Visualization ===
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

for atom in atoms:
    symbol, x, y, z = atom
    color = atom_colors.get(symbol, "gray")   # default gray
    size = atom_sizes.get(symbol, 150)        # default size
    ax.scatter(x, y, z, color=color, s=size, label=symbol)

# Optional: enable atom labels
# for atom in atoms:
#     symbol, x, y, z = atom
#     ax.text(x, y, z, symbol, fontsize=8, color="black")

# === Axis settings ===
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title("Molecule Visualization (Atoms by Type & Size)")

# Prevent duplicate legend entries
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
ax.legend(unique.values(), unique.keys())

plt.tight_layout()
plt.show()
