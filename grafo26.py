#!/usr/bin/env python3
import os
import shutil
import string

# ============================================================
# CONFIGURACIÓN GENERAL
# ============================================================

BASE = "FINAL_HL_C6H6"
RXNET = os.path.join(BASE, "RXNet.cg")
RXNET_BL = os.path.join(BASE, "RXNet.barrless")

OUTDIR = "OUTPUT"
PILGRIM_DIR = "UDATA"
CHEM_DIR = "chem_C6H6"
CHEM_FILE = os.path.join(CHEM_DIR, "pif.chem")

TS_BARRIERLESS_OFFSET = 800

# ============================================================
# ENTRADAS DEL USUARIO
# ============================================================

fmin  = "FINAL_HL_C6H6/MINinfo_zpecorr"
fprod = "FINAL_HL_C6H6/PRODinfo_zpecorr"
fts   = "FINAL_HL_C6H6/TSinfo_zpecorr"

start_path = ["MIN1"]
max_length = 3 #Cuantos min hay entre min1 y prx. Ejemplo: MIN1(0.0) --TS4(60.1)--> MIN13(56.0) --TS804(93.7)--> PR50(93.7) Aquí hay 1 min intermedio, entonces con max_length de 1 apareceria.
E_cutoff = 142.7 #kcal/mol

# ============================================================
# COLORES ANSI
# ============================================================

RESET   = "\033[0m"
BOLD    = "\033[1m"
RED     = "\033[91m"
BLUE    = "\033[94m"
GREEN   = "\033[92m"
MAGENTA = "\033[95m"
GRAY    = "\033[90m"

# ============================================================
# CABECERA pif.chem
# ============================================================

CHEM_HEADER = """#-----------------------------------------------#
#  Pilgrim pif.chem
#-----------------------------------------------#
"""

# ============================================================
# FUNCIÓN PARA GENERAR ETIQUETAS a, b, ..., z, a1, ..., z1, a2...
# ============================================================

def reaction_label(i):
    letters = string.ascii_lowercase
    letter = letters[i % 26]
    cycle = i // 26
    if cycle == 0:
        return letter
    else:
        return f"{letter}{cycle}"

# ============================================================
# LECTURA DE ENERGÍAS
# ============================================================

def read_energies(fname, prefix):
    energies = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                idx = int(parts[0])
                val = float(parts[1])
            except:
                continue
            energies[f"{prefix}{idx}"] = val
    return energies

def load_all_energies():
    e = {}
    e.update(read_energies(fmin,  "MIN"))
    e.update(read_energies(fprod, "PR"))
    e.update(read_energies(fts,   "TS"))
    return e

# ============================================================
# PARSER RXNet
# ============================================================

def parse_rxnet(fname, barrierless=False):
    reactions = []

    if not os.path.isfile(fname):
        return reactions

    with open(fname) as f:
        for line in f:
            l = line.split()
            if not l or not l[0].isdigit() or "DIS" in line:
                continue

            tsnum = int(l[0])
            if barrierless:
                tsnum += TS_BARRIERLESS_OFFSET

            ts_label = f"TS{tsnum}"
            rmin = f"MIN{int(l[3])}"

            # ← NUEVO: detectar tipo de flecha
            bidirectional = (l[4] == "<--->")

            if "PR" in l[5]:
                prod = f"PR{int(l[5][2:-1])}"
            else:
                prod = f"MIN{int(l[6])}"

            # Dirección directa — siempre
            reactions.append((rmin, ts_label, prod, barrierless))

            # ← NUEVO: dirección inversa si es <---> y prod es MIN
            if bidirectional and not prod.startswith("PR"):
                reactions.append((prod, ts_label, rmin, barrierless))

    return reactions

def build_reactions():
    return parse_rxnet(RXNET, False) + parse_rxnet(RXNET_BL, True)

# ============================================================
# FILTRO ENERGÉTICO
# ============================================================

def path_is_valid(path, edge2ts, energies):

    if len(path) - 2 > max_length:   # descuenta MIN1 y el PR final
        return False

    for sp in path:
        if energies.get(sp, 0.0) > E_cutoff:
            return False

    for i in range(len(path) - 1):
        ts = edge2ts[(path[i], path[i+1])]
        if energies.get(ts, 0.0) > E_cutoff:
            return False

    return True

# ============================================================
# BÚSQUEDA DE CAMINOS
# ============================================================

def find_paths(reactions, energies):

    adj = {}
    edge2ts = {}
    edge2nobar = {}

    unique_reactions = []
    seen_reactions = set()

    for r, ts, p, nobar in reactions:
        key = (r, ts, p, nobar)
        if key not in seen_reactions:
            seen_reactions.add(key)
            unique_reactions.append((r, ts, p, nobar))

    for r, ts, p, nobar in unique_reactions:
        e_ts = energies.get(ts, 0.0)
        if (r, p) not in edge2ts:
            adj.setdefault(r, []).append(p)
            edge2ts[(r, p)] = ts
            edge2nobar[(r, p)] = nobar
        else:
        # Conservar el TS con menor energía
            if e_ts < energies.get(edge2ts[(r, p)], 0.0):
                edge2ts[(r, p)] = ts
                edge2nobar[(r, p)] = nobar
        if nobar:
            energies[ts] = energies.get(p, 0.0)
    valid_paths = []

    def dfs(path):
        last = path[-1]

        if last.startswith("PR"):
            if path_is_valid(path, edge2ts, energies):
                valid_paths.append(path)
            return

        if len(path) >= max_length + 2:   # MIN1 + intermedios + PR
            return

        for nxt in adj.get(last, []):
            if nxt in path:
                continue
            dfs(path + [nxt])

    dfs(start_path)

    unique_paths = []
    seen_paths = set()

    for p in valid_paths:
        t = tuple(p)
        if t not in seen_paths:
            seen_paths.add(t)
            unique_paths.append(p)

    return unique_paths, edge2ts, edge2nobar

# ============================================================
# IMPRESIÓN (NO SE TOCA)
# ============================================================

def print_paths(paths, reactions, energies, edge2ts, edge2nobar):

    paths = sorted(paths, key=lambda p: int(p[-1][2:]))

    print("\n" + BOLD + "===== CAMINOS REACCIONALES VÁLIDOS =====" + RESET + "\n")

    total_paths = len(paths)
    paths_with_nobar = 0

    all_reactions = []
    barrierless_count = 0
    barrier_count = 0
    unique_barrier = set()

    for i, path in enumerate(paths, 1):

        has_nobar = False
        pieces = []

        for j, sp in enumerate(path):

            if sp.startswith("MIN"):
                color_sp = GREEN
            elif sp.startswith("PR"):
                color_sp = MAGENTA
            else:
                color_sp = RESET

            pieces.append(f"{color_sp}{sp}({energies.get(sp,0):.1f}){RESET}")

            if j < len(path)-1:
                ts = edge2ts[(sp, path[j+1])]
                nobar = edge2nobar[(sp, path[j+1])]

                if nobar:
                    color_ts = RED
                    barrierless_count += 1
                    has_nobar = True
                else:
                    color_ts = BLUE
                    barrier_count += 1
                    unique_barrier.add((sp, ts, path[j+1]))

                all_reactions.append((sp, ts, path[j+1], nobar))
                pieces.append(f"--{color_ts}{ts}({energies.get(ts,0):.1f}){RESET}-->")

        if has_nobar:
            paths_with_nobar += 1

        print(GRAY + "-"*70 + RESET)
        print(BOLD + f"Camino {i}  →  {path[-1]}" + RESET)
        print(" ".join(pieces))
        print()

    print(GRAY + "="*70 + RESET)
    print(BOLD + "===== RESUMEN GLOBAL =====" + RESET)
    print(f"  Caminos totales : {total_paths}  ({paths_with_nobar} con alguna reacción sin barrera)")
    print(f"  Reacciones      : {len(all_reactions)}  ({barrier_count} con barrera | {len(unique_barrier)} únicas | {barrierless_count} sin barrera)")

    # ---- Caminos por producto ----
    from collections import Counter
    prod_count = Counter(p[-1] for p in paths)
    print()
    print(GRAY + "-"*70 + RESET)
    print("  Caminos por producto:")
    for prod, count in sorted(prod_count.items(), key=lambda x: int(x[0][2:])):
        print(f"    {MAGENTA}{prod}{RESET}  →  {count} camino{'s' if count > 1 else ''}")

    # ---- Camino(s) de mínima energía máxima ----
    def max_energy_path(path):
        e_nodes = [energies.get(sp, 0.0) for sp in path]
        e_ts    = [energies.get(edge2ts[(path[i], path[i+1])], 0.0) for i in range(len(path)-1)]
        return max(e_nodes + e_ts)

    min_emax = min(max_energy_path(p) for p in paths)
    best_paths = [p for p in paths if max_energy_path(p) == min_emax]

    print()
    print(GRAY + "-"*70 + RESET)
    print(f"  Camino(s) de mínima energía máxima (cuello de botella más bajo):")
    for bp in best_paths:
        pieces = []
        for j, sp in enumerate(bp):
            color_sp = GREEN if sp.startswith("MIN") else MAGENTA
            pieces.append(f"{color_sp}{sp}({energies.get(sp,0):.1f}){RESET}")
            if j < len(bp)-1:
                ts = edge2ts[(sp, bp[j+1])]
                nobar = edge2nobar[(sp, bp[j+1])]
                color_ts = RED if nobar else BLUE
                pieces.append(f"--{color_ts}{ts}({energies.get(ts,0):.1f}){RESET}-->")
        print(f"    {' '.join(pieces)}")
        print(f"    [ Energía máxima del camino: {min_emax:.1f} kcal/mol ]")

    # ---- Criterio de cutoff ----
    margen = E_cutoff - min_emax
    print()
    print(GRAY + "-"*70 + RESET)
    print(f"  Criterio de corte aplicado:")
    print(f"    E_cutoff                        = {E_cutoff:.1f} kcal/mol")
    print(f"    Energía máxima del camino menos energético = {min_emax:.1f} kcal/mol")
    print(f"    Margen disponible               = {E_cutoff:.1f} - {min_emax:.1f} = {margen:.1f} kcal/mol")
    print(f"    → En combustión los caminos de menor barrera dominan cinéticamente.")
    print(f"      Todos los caminos mostrados tienen su paso limitante por debajo")
    print(f"      del cutoff energético definido ({E_cutoff:.1f} kcal/mol).")
    print(GRAY + "="*70 + RESET + "\n")
    return paths

# ============================================================
# ESCRITURA pif.chem (con barrera + sin barrera juntos)
# ============================================================

def write_pif_chem(valid_paths, edge2ts, edge2nobar):

    os.makedirs(CHEM_DIR, exist_ok=True)

    valid_paths = sorted(valid_paths, key=lambda p: int(p[-1][2:]))

    written_barrier = set()
    written_nobar = set()

    with open(CHEM_FILE, "w") as f:
        f.write(CHEM_HEADER)

        # -------- CON BARRERA (numeradas) --------
        idx = 1
        for path in valid_paths:
            for i in range(len(path)-1):
                r = path[i]
                p = path[i+1]
                ts = edge2ts[(r, p)]
                nobar = edge2nobar[(r, p)]

                if nobar:
                    continue

                key = (r, ts, p)
                if key in written_barrier:
                    continue

                if p.startswith("PR"):
                    p_written = f"{p}_a + {p}_b"
                else:
                    p_written = p

                f.write(f"{idx} : {r} --> {ts} --> {p_written}\n")
                written_barrier.add(key)
                idx += 1

        # -------- SIN BARRERA (letras + TS_dummy) --------
        f.write("\n")
        letter_idx = 0
        for path in valid_paths:
            for i in range(len(path)-1):
                r = path[i]
                p = path[i+1]
                nobar = edge2nobar[(r, p)]

                if not nobar:
                    continue

                key = (r, p)
                if key in written_nobar:
                    continue

                label = reaction_label(letter_idx)
                ts_dummy = f"TS_dummy_{label}"

                if p.startswith("PR"):
                    p_written = f"{p}_a + {p}_b"
                else:
                    p_written = p

                f.write(f"{label} : {r} --> {ts_dummy} --> {p_written}\n")
                written_nobar.add(key)
                letter_idx += 1

        total_unique_nobar = len(written_nobar)
        f.write(f"\n# Total barrierless reactions (unique): {total_unique_nobar}\n")

    print(f"Reacciones sin barrera únicas: {total_unique_nobar}")
# ============================================================
# COPIA DE ESTRUCTURAS
# ============================================================

def collect_structures(valid_paths, edge2ts):

    os.makedirs(PILGRIM_DIR, exist_ok=True)

    species = set()

    for path in valid_paths:
        for sp in path:
            species.add(sp)
        for i in range(len(path)-1):
            species.add(edge2ts[(path[i], path[i+1])])

    for sp in species:
        if sp.startswith("PR"):
            files = [f"{sp}_a.out", f"{sp}_b.out"]
        else:
            files = [f"{sp}.out"]

        for f in files:
            src = os.path.join(OUTDIR, f)
            if os.path.isfile(src):
                shutil.copy(src, os.path.join(PILGRIM_DIR, f))

# ============================================================
# MAIN
# ============================================================

def main():
    reactions = build_reactions()
    energies = load_all_energies()

    valid_paths, edge2ts, edge2nobar = find_paths(reactions, energies)
    ordered_paths = print_paths(valid_paths, reactions, energies, edge2ts, edge2nobar)
    write_pif_chem(ordered_paths, edge2ts, edge2nobar)
    collect_structures(ordered_paths, edge2ts)
    print("Proceso completado correctamente ✔\n")

if __name__ == "__main__":
    main()

