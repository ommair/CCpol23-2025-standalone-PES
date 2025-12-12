import numpy as np
import pandas as pd

# === Atom coordinates in bohr from your data ===
coords = {
    "O":  [0.00000000000000, 0.0000000, 0.1255334885],
    "H1": [-1.45365196228170, 0.0000000, -0.9961538357],
    "H2": [ 1.45365196228170, 0.0000000, -0.9961538357],
    "D1a": [0.00000000000000, -0.2067213, -0.2462589474],
    "D1b": [0.00000000000000,  0.2067213, -0.2462589474],
    "N1a": [ 0.90,  0.053,  0.64],
    "N1b": [-0.90,  0.053,  0.64],
    "N1c": [ 0.90, -0.053,  0.64],
    "N1d": [-0.90, -0.053,  0.64],
    "N2a": [ 1.29,  0.14, -0.91],
    "N2b": [-1.29,  0.14, -0.91],
    "N2c": [ 1.29, -0.14, -0.91],
    "N2d": [-1.29, -0.14, -0.91],
    "N3a": [ 0.56,  0.24, -0.48],
    "N3b": [-0.56,  0.24, -0.48],
    "N3c": [ 0.56, -0.24, -0.48],
    "N3d": [-0.56, -0.24, -0.48],
    "N4a": [ 0.91,  0.21, -0.68],
    "N4b": [-0.91,  0.21, -0.68],
    "N4c": [ 0.91, -0.21, -0.68],
    "N4d": [-0.91, -0.21, -0.68],
    "N5a": [ 1.48,  0.26, -0.62],
    "N5b": [-1.48,  0.26, -0.62],
    "N5c": [ 1.48, -0.26, -0.62],
    "N5d": [-1.48, -0.26, -0.62],
}


# === Define local frame: O, H1, H2
a = np.array(coords["O"])
b = np.array(coords["H1"])
c = np.array(coords["H2"])

ba = a - b
ca = a - c
dd = np.cross(ba, ca)  # right-handed orthonormal basis

# === Solve for coefficients ===
def get_c123(rD):
    rhs = np.array(rD) - a
    M = np.column_stack((ba, ca, dd))
    return np.linalg.solve(M, rhs)

# === Reconstruct from coefficients
rows = []
max_error = 0
for label, rD in coords.items():
    if label in ["O", "H1", "H2"]:
        continue  # Skip base atoms
    c1, c2, c3 = get_c123(rD)
    rD_calc = a + c1 * ba + c2 * ca + c3 * dd
    error = np.linalg.norm(np.array(rD) - rD_calc)
    max_error = max(max_error, error)
    rows.append({
        "site": label,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "err_x": rD[0] - rD_calc[0],
        "err_y": rD[1] - rD_calc[1],
        "err_z": rD[2] - rD_calc[2],
        "rmsd": error
    })

df = pd.DataFrame(rows)
print(df[["site", "c1", "c2", "c3", "rmsd"]].to_string(index=False))
print(f"\nMax RMSD = {max_error:.3e} ")

