# Control rod banked height calculation
# Reuven Z. Lazarus, 2013
# rlazarus@alumni.reed.edu

import subprocess

import numpy as np

import rodcal

POWERS = [0, 20, 40, 50, 60, 80, 100, 120, 140, 150, 160, 180, 200, 220, 230]
CORE_EXCESSES = [0.05 * i for i in xrange(0, 41)]  # 0.00, 0.05, ..., 2.00

def rod_height_vs_core_excess(rod_banks):
    """
    Calculate the cubic of best fit for the banked rod heights vs. core excess.

    Args:
        rod_banks: Data collected from the reactor. A list of tuples of
            (power in kW, banked height in %, core excess in $).

    Returns:
        The cubic of best fit, as a poly1d instance.
    """
    x = [t[2] for t in rod_banks]
    y = [t[1] for t in rod_banks]
    return np.poly1d(np.polyfit(x, y, 3))

def core_excess_vs_power(rod_banks):
    """
    Calculate the cubic of best fit for the core excess vs. power.

    Args:
        rod_banks: Data collected from the reactor. A list of tuples of
            (power in kW, banked height in %, core excess in $).

    Returns:
        The cubic of best fit, as a poly1d instance.
    """
    x = [t[0] for t in rod_banks]
    y = [t[2] for t in rod_banks]
    return np.poly1d(np.polyfit(x, y, 3))

def rod_height(core_excess, power, hvsrho, rhovsp):
    """
    Calculate the target rod height for a given core excess and power.

    Args:
        core_excess: The five-watt core excess, in dollars.
        power: The target power.
        hvsrho: The output of rod_height_vs_core_excess.
        rhovsp: The output of core_excess_vs_power.

    Returns:
        The rod height in %.
    """
    reactivity = rhovsp(power) + core_excess - rhovsp(0)
    return hvsrho(reactivity)

def table(hvsrho, rhovsp):
    """
    Generate the table of banked rod heights for all core excesses and powers.

    Args:
        hvsrho: The output of rod_height_vs_core_excess.
        rhovsp: The output of core_excess_vs_power.

    Returns:
        [[rod height for each power] for each core excess]
    """
    return [[rod_height(core_excess, power, hvsrho, rhovsp)
             for power in POWERS]
            for core_excess in CORE_EXCESSES]

def tex(xenon_free, data):
    """
    Generate and typeset a printable table of estimated banked rod heights for
    the back of the logbook.

    Args:
        xenon_free: The experimentally-determined 5 W core excess without xenon.
        data: The output of table().

    Returns:
        The filename of the final typeset PDF.
    """
    # We'll want to bold the rows near the xenon-free core excess
    with open("template-banktable.tex") as f:
        tex = f.read()

    replacements = {}
    replacements["numpowers"] = len(POWERS)
    replacements["powers"] = " & ".join(r"\textbf{" + str(x) + "}"
                                        for x in POWERS)
    replacements["tablebody"] = ""

    for i, core_excess in enumerate(CORE_EXCESSES):
        bold = (abs(core_excess - xenon_free) <= 0.05)

        if bold:
            replacements["tablebody"] += r"\textbf{"
        replacements["tablebody"] += r"\${:.2f}".format(core_excess)
        if bold:
            replacements["tablebody"] += "}"

        for height in data[i]:
            replacements["tablebody"] += " & "

            if height < 100.5:
                if bold:
                    replacements["tablebody"] += r"\textbf{"
                replacements["tablebody"] += str(int(round(height)))
                if bold:
                    replacements["tablebody"] += "}"

        replacements["tablebody"] += r"\\"

    tex = tex % replacements

    filename = rodcal.available_filename("banktable")
    with open(filename, 'w') as f:
        f.write(tex)

    subprocess.call(["pdflatex", filename, "-interaction", "batchmode"])

    return filename.replace("tex", "pdf")

def main():
    # Data entry
    print "For each power, enter the banked rod height and core excess from"
    print "Page 4 of SOP 34A (Control Rod Calibration Form)."
    print

    rod_banks = []
    for power in [0.05, 1, 25, 50, 75, 100, 125, 150, 175, 200, 225, 230]:
        print "          Power: {} kW".format(power)
        height = float(raw_input(" Rod height (%): "))
        cxs = float(raw_input("Core excess ($): "))
        rod_banks.append((power, height, cxs))
        print

    hvsrho = rod_height_vs_core_excess(rod_banks)
    rhovsp = core_excess_vs_power(rod_banks)
    xenon_free = rod_banks[0][2]

    filename = tex(xenon_free, table(hvsrho, rhovsp))
    print "Filename: " + filename
    rodcal.open_file(filename)

if __name__ == "__main__":
    main()
