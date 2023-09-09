import re
import os

class GamessOut:
    def __init__(self):
        self.success = False
        self.error_message = ""
        self.total_energy = None
        self.coordinates = []

    def __str__(self):
        return "Success:{}".format(self.success)


def gparse(gamout, parse_type="default"):
    r = GamessOut()
    out_str = open(gamout, "r").read()

    if os.name == "nt":
        r.success = True
    elif os.name == "posix":
        # exited gracefully?
        if not out_str.endswith("gracefully.\n"):
            r.error_message = out_str[-1000:]
            return r
        else:
            r.success = True
    else:
        r.success = True

    # NOTE: Above check appears to be version-dependent.
    #   Temporary fix: assume it worked until proven otherwise:

    if parse_type == "default":
        r = default_parse(out_str, r)

    return r


def default_parse(out_str, r):
    eng_re = re.compile('TOTAL ENERGY =(.*)\n')
    pcm_eng_re = re.compile('RESULTS OF PCM CALCULATION(.*?)A\.U\.\n\n', re.DOTALL)
    coord_re = re.compile('COORDINATES OF ALL ATOMS ARE (.*?)------------\n(.*?)\n\n', re.DOTALL)
    pop_re = re.compile('TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS(.*?)\n\n',re.DOTALL)
    es_moment_re = re.compile('DEBYE\)(.*?)\n \.\.', re.DOTALL)
    num_electron_re = re.compile('NUMBER OF ELECTRONS(.*?)\n')
    mo_re = re.compile('MOLECULAR ORBITALS(.*?)\n\n     ------', re.DOTALL)
    eigen_re = re.compile('EIGENVECTORS(.*?)\n\n     ------', re.DOTALL)
    hessian_re = re.compile('THE HARMONIC ZERO POINT ENERGY IS(.*?)KCAL/MOL', re.DOTALL)
    ir_re = re.compile('SYMMETRY  RED\. MASS  IR INTENS\.\n(.*?)\n\n     ------', re.DOTALL)
    stationary_point_re = re.compile('THIS IS NOT A STATIONARY POINT ON THE MOLECULAR PES')
    nsearch_re = re.compile('    NSERCH=(.*)\n')
    tddft_re = re.compile('SUMMARY OF TDDFT RESULTS\n\n(.*?)\n TRANSITION', re.DOTALL)
    nmr_re = re.compile('GIAO CHEMICAL SHIELDING TENSOR \(PPM\):\n(.*?)DONE WITH NMR SHIELDINGS', re.DOTALL)

    # Total Energy, this only match in gas phase calculations
    r.total_energy = None
    for m in eng_re.finditer(out_str):
        r.total_energy = float(m.group(1).strip())
    # parse PCM result if eng_re didn't match
    m = pcm_eng_re.search(out_str)
    if m is not None:
        for l in m.group().split("\n"):
            if l.endswith("A.U."):
                if l.startswith(" FREE ENERGY IN SOLVENT"):
                    r.free_energy = float(l[-20:-5].strip())
                elif l.startswith(" INTERNAL ENERGY IN SOLVENT"):
                    r.internal_energy = float(l[-20:-5].strip())
                elif l.startswith(" DELTA INTERNAL ENERGY"):
                    r.delta_internal_energy = float(l[-20:-5].strip())
                elif l.startswith(" ELECTROSTATIC INTERACTION"):
                    r.electrostatic_interaction = float(l[-20:-5].strip())
                elif l.startswith(" PIEROTTI CAVITATION ENERGY"):
                    r.pierotti_cavitation_energy = float(l[-20:-5].strip())
                elif l.startswith(" DISPERSION FREE ENERGY"):
                    r.dispersion_free_energy = float(l[-20:-5].strip())
                elif l.startswith(" REPULSION FREE ENERGY"):
                    r.repulsion_free_energy = float(l[-20:-5].strip())
                elif l.startswith(" TOTAL INTERACTION"):
                    r.total_interacion = float(l[-20:-5].strip())
                elif l.startswith(" TOTAL FREE ENERGY"):
                    r.total_energy = float(l[-20:-5].strip())
                else:
                    pass

    # NSERCH (optimization enegy transition)
    r.nsearches = []
    for m in nsearch_re.finditer(out_str):
        l = m.group(1).strip()
        r.nsearches.append(float(l.split()[-1]))

    # Coordinates
    for m in coord_re.finditer(out_str):
        r.coordinates = []
        for l in m.group(2).split("\n"):
            coord = l.split()
            cds = [float(coord[2]), float(coord[3]), float(coord[4])]
            r.coordinates.append(cds)

    # MULLIKEN and LOWDIN charge
    for m in pop_re.finditer(out_str):
        r.mulliken_populations = []
        r.lowdin_populations = []
        r.mulliken_charges = []
        r.lowdin_charges = []
        for l in m.group(1).split("\n"):
            ls = l.split()
            if len(ls) == 6:
                r.mulliken_populations.append(float(ls[2]))
                r.mulliken_charges.append(float(ls[3]))
                r.lowdin_populations.append(float(ls[4]))
                r.lowdin_charges.append(float(ls[5]))

    # Dipole morment
    for m in es_moment_re.finditer(out_str):
        r.dipole_moment = []
        for l in m.group(1).split("\n"):
            ls = l.split()
            if len(ls) == 4:
                for v in ls:
                    r.dipole_moment.append(float(v))

    # Num of electrons
    m = num_electron_re.search(out_str)
    if m is not None:
        num_elec = m.group().split("=")[1][:-1]
        num_elec = int(num_elec)

    # MO
    m = eigen_re.search(out_str)
    if m is not None:
        mo_energies = []
        for l in m.group(1).split("\n"):
            # startwith many spaces and line contains float value
            if l.startswith("                  ") and l.find(".") > 0:
                ls = [float(v) for v in l.split()]
                mo_energies += ls

    r.orbital_energies = mo_energies
    r.nHOMO = r.orbital_energies[int(num_elec/2)-2]
    r.HOMO = r.orbital_energies[int(num_elec/2)-1]
    r.LUMO = r.orbital_energies[int(num_elec/2)]
    r.nLUMO = r.orbital_energies[int(num_elec/2)+1]

    # Optimization output file contains "MOLECULAR ORBITALS" section.
    m = mo_re.search(out_str)
    if m is not None:
        mo_energies = []
        for l in m.group(1).split("\n"):
            # startwith many spaces and line contains float value
            if l.startswith("                  ") and l.find(".") > 0:
                ls = [float(v) for v in l.split()]
                mo_energies += ls

    # HESSIAN
    r.is_stationary_point = None
    r.ZPE = None
    m = hessian_re.search(out_str)
    if m is not None:
        sp = stationary_point_re.search(out_str)
        if sp is not None:
            r.is_stationary_point = False
        else:
            r.ZPE =float(m.group(1).split("\n")[-1])

    ir_spectra = []
    m = ir_re.search(out_str)
    if m is not None:
        for l in m.group(1).split("\n"):
            ls = l.split()
            ir_spectra.append((float(ls[1]), float(ls[4])))
        r.ir_spectra = ir_spectra


    # TDDFT
    uv_spectra = []
    m = tddft_re.search(out_str)
    if m is not None:
        for l in m.group(1).split("\n")[2:]:
            ls = l.split()
            if len(ls) == 8:
                uv_spectra.append((float(ls[3]), float(ls[7])))
        r.uv_spectra = uv_spectra

    # NMR
    isotropic_shielding = []
    m = nmr_re.search(out_str)
    if m is not None:
        for l in m.group(1).split("\n"):
            ls = l.split()
            if len(ls) == 1:
                try:
                    n = float(ls[0])
                    isotropic_shielding.append(n)
                except:
                    pass
        r.isotropic_shielding = isotropic_shielding

    return r
