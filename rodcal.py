# Control rod calibration
# Reuven Z. Lazarus, 2013
# rlazarus@alumni.reed.edu

import math
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt

# Constants

beta = 0.0065
beta_eff = 0.0075
l_p = 0.000108  # seconds

# values for the six delayed neutron groups - ref. Lamarsh, Table 3.5
t_half = [55.72, 22.72, 6.22, 2.30, 0.610, 0.230]
lambda_i = [math.log(2)/x for x in t_half]
beta_i = [0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]


def reactivity_calc(bottom, pulls):
    """
    Perform the actual rod worth calculations for a single control rod.

    Args:
        bottom: The indicated rod height for this rod at the lowest point during
            the calibration.
        pulls: A list of (height, period): the height at the top of the pull and
            the stable period observed following it.

    Returns:
        A 4-tuple: (integral rod worth polynomial, integral rod worth data,
                    reactivity addition rate polynomial, addition rate data).
        Polynomials are instances of np.poly1d. Data are lists of tuples.
    """
    # Calculate the reactivity insertion for each rod pull
    reactivity = []  # (%height start, %height end, $)
    for (i, (end, period)) in enumerate(pulls):
        if i == 0:
            start = bottom
        else:
            start = pulls[i - 1][0]

        reactivity_dkk = (l_p/period +
                          sum(beta / (1.0 + lambda_ * period) for beta, lambda_
                              in zip(beta_i, lambda_i)))
        reactivity_dollar = reactivity_dkk / beta_eff
        reactivity.append((start, end, reactivity_dollar))

    # Calculate the integral rod worth after each pull, and the average
    # reactivity addition rate during each pull.
    integral = []  # (%height, $)
    addition_rate = []  # (%height, $/%)

    integral = [(reactivity[i][0], sum(reactivity[j][2] for j in
                                       xrange(i, len(reactivity))))
                for i in xrange(len(reactivity))]

    addition_rate = [((reactivity[i][0] + reactivity[i][1]) / 2.0,
                      reactivity[i][2]/(reactivity[i][1] - reactivity[i][0]))
                     for i in xrange(len(reactivity))]

    # Clamp the integral worth to zero when fully withdrawn
    integral.append((100.0, 0.0))
    # Clamp the addition rate to zero at the top and bottom
    addition_rate.insert(0, (0.0, 0.0))
    addition_rate.append((100.0, 0.0))

    # Fit the integral worth to a cubic
    x, y = zip(*integral)
    integral_fit = np.poly1d(np.polyfit(x, y, 3))

    # Fit the addition rate to a cubic
    x, y = zip(*addition_rate)
    addition_rate_fit = np.poly1d(np.polyfit(x, y, 3))

    return (integral_fit, integral, addition_rate_fit, addition_rate)

def plot(rod, integral_fit, integral, addition_rate_fit, addition_rate):
    """
    Generate graphs from the integral fit and reactivity addition rate.

    Args:
        rod: "safe", "shim", or "reg"
        integral_fit: A np.poly1d polynomial.
        integral: A list of tuples: the raw data from which the fit was computed.
        addition_rate_fit: A np.poly1d polynomial.
        addition_rate: Data for the fit.
    """
    r = np.arange(0,100,0.1)

    plt.figure()
    # Plot the fit in a blue line
    plt.plot(r, np.poly1d(integral_fit)(r), "b-")
    # Plot the data in black circles
    x, y = zip(*integral)
    plt.plot(x, y, "ko")

    plt.xlabel("Rod Height (%)")
    plt.ylabel("Integral Worth ($)")
    plt.tight_layout(pad=0)
    plt.savefig("{}-integral.png".format(rod))

    plt.figure()
    plt.plot(r, np.poly1d(addition_rate_fit)(r), "b-")

    x, y = zip(*addition_rate)
    plt.plot(x, y, "ko")

    plt.xlabel("Rod Height (%)")
    plt.ylabel("Insertion Rate ($/%)")
    plt.tight_layout(pad=0)
    plt.savefig("{}-rate.png".format(rod))

def tech_specs(safe_int, safe_add, shim_int, shim_add, reg_int, reg_add,
               safe_critical, shim_critical, withdrawal_times):
    """
    Check for Tech Specs compliance and do whatever other reporting and
    summarizing is necessary.

    Args:
        safe_int: The safe rod's integral rod worth fit, as a numpy poly1d.
        safe_add: The safe rod's reactivity addition rate fit, as a poly1d.
        shim_int: The shim rod's integral rod worth fit.
        shim_add: The shim rod's reactivity addition rate fit.
        reg_int: The reg rod's integral rod worth fit.
        reg_add: The reg rod's reactivity addition rate fit.
        safe_critical: The lowest critical safe rod height.
        shim_critical: The lowest critical shim rod height.
        withdrawal_times: List containing the time it takes to pull the [safe,
            shim, reg] rods all the way out, in seconds.

    Returns:
        A dict of the values to go into the report.
    """
    result = {}

    # TeXify the polynomials
    for key, poly in [("safeintpoly", safe_int),
                      ("safeaddpoly", safe_add),
                      ("shimintpoly", shim_int),
                      ("shimaddpoly", shim_add),
                      ("regintpoly", reg_int),
                      ("regaddpoly", reg_add)]:
        patterns = ["{:.6} \, x^3", "{:.6} \, x^2", "{:.6} \, x", "{:.6} "]
        strings = [p.format(poly.c[i]) for i,p in enumerate(patterns)]
        for i in xrange(len(strings)):
            strings[i] = strings[i].format(poly.c[i])
            # Render the scientific notation nicely
            strings[i] = re.sub("(.+)e([+-])0*([^ ]+) ",
                                r"(\1 \\times 10^{\2\3})",
                                strings[i])
        result[key] = "$" + " + ".join(strings) + "$"

    # Find the total rod worths
    result["safeworth"] = safe_int(0.0)
    result["shimworth"] = shim_int(0.0)
    result["regworth"] = reg_int(0.0)
    result["totalworth"] = (result["safeworth"] + result["shimworth"] +
                            result["regworth"])

    # Find the highest reactivity addition rate ($/%)
    # Find the extremes of the addition-rate curve; take the one in range
    clamp = lambda x: 0 <= x <= 100
    result["safemaxdpp"] = safe_add(filter(clamp, safe_add.deriv().r)[0]) * 100
    result["shimmaxdpp"] = shim_add(filter(clamp, shim_add.deriv().r)[0]) * 100
    result["regmaxdpp"] = reg_add(filter(clamp, reg_add.deriv().r)[0]) * 100

    # Convert to cent/s using the withdrawal times
    result["safemaxdps"] = result["safemaxdpp"] * 100 / withdrawal_times[0]
    result["shimmaxdps"] = result["shimmaxdpp"] * 100 / withdrawal_times[1]
    result["regmaxdps"] = result["regmaxdpp"] * 100 / withdrawal_times[2]

    # Check against TS maximum 12c/sec
    result["safedpsok"] = result["safemaxdps"] < 12
    result["shimdpsok"] = result["shimmaxdps"] < 12
    result["regdpsok"] = result["regmaxdps"] < 12

    # Find the CXS independently using the Safe and Shim rods
    result["safecxsht"] = safe_critical
    result["shimcxsht"] = shim_critical

    result["safecxs"] = safe_int(safe_critical)
    result["shimcxs"] = shim_int(shim_critical)

    # Check against TS maximum $3.00
    result["cxsok"] = result["safecxs"] < 3 and result["shimcxs"] < 3

    # Find the shutdown margin
    result["safesdm"] = result["totalworth"] - result["safecxs"]
    result["shimsdm"] = result["totalworth"] - result["shimcxs"]

    # Check against TS minimum $1.00 to shut down
    result["sdmok"] = result["safesdm"] > 1 and result["shimsdm"] > 1

    # Find the one-stuck-rod shutdown margin
    if result["safeworth"] > result["shimworth"]:
        result["mostrxvrod"] = "Safe Rod"
        stuckrod = result["safeworth"]
    else:
        result["mostrxvrod"] = "Shim Rod"
        stuckrod = result["shimworth"]

    result["safeosr"] = result["totalworth"] - result["safecxs"] - stuckrod
    result["shimosr"] = result["totalworth"] - result["shimcxs"] - stuckrod

    # Check against TS minimum $0.50
    result["osrok"] = result["safeosr"] > 0.5 and result["shimosr"] > 0.5

    return result

def tex_report(result):
    """
    Generate and typeset a printable report of the calibration, including rod
    worth summary and Tech Specs check. The inputs aren't reproduced; this is
    meant to be attached to the worksheet, not to replace it.

    Args:
        result: The results of the calibration; the dict returned by
            tech_specs().
    Returns:
        The filename of the final typeset PDF.
    """
    # Format values for printing
    # Dollars
    for key in ("safeworth", "shimworth", "regworth", "totalworth", "safecxs",
                "shimcxs", "safesdm", "shimsdm", "safeosr", "shimosr"):
        result[key] = "\\${:.2f}".format(result[key])

    # Rod heights
    for key in "safecxsht", "shimcxsht":
        result[key] = "{:.1f}\\%".format(result[key])

    # Insertion rates
    for key in ("safemaxdpp", "shimmaxdpp", "regmaxdpp"):
        result[key] = "{:.2f} \\textcentoldstyle/\\%".format(result[key])
    for key in ("safemaxdps", "shimmaxdps", "regmaxdps"):
        result[key] = "{:.2f} \\textcentoldstyle/sec".format(result[key])

    # Booleans
    for key in "safedpsok", "shimdpsok", "regdpsok", "cxsok", "sdmok", "osrok":
        if result[key]:
            result[key] = "OK"
        else:
            result[key] = "PROBLEM"

    # And "mostrxvrod" is a string, so we leave it as is.
    with open("template-report.tex") as f:
        tex = f.read()

    tex = tex % result

    filename = available_filename("report")
    with open(filename, 'w') as f:
        f.write(tex)

    subprocess.call(["pdflatex", filename, "-interaction", "batchmode"])

    return filename.replace("tex", "pdf")

def tabular_data(cubics):
    """
    Transforms a list of three cubics into lists of height-worth tuples for the
    tables in the back of the logbook. Safe and Shim tables start from 48%, Reg
    table from 0%.

    Args:
        cubics: [safe, shim, reg], where each is the integral rod worth fit, as
            a numpy poly1d instance.

    Returns:
        [safe, shim, reg], where each is a list of (%, $) increasing by 0.1%
            from 0.0% (for the reg) or 48.0% (for the Safe and Shim) to 100.0%.
    """
    result = []
    for n, p in enumerate(cubics):
        range = np.arange(0.0 if n == 2 else 48.0, 100.0, 0.1)
        table = zip(range, p(range))
        result.append(table)
    return result

def tex_tables(tabular_data, safe_worth, shim_worth, reg_worth):
    """
    Generate and typeset a printable set of rod worth tables for the back of
    the logbook.

    Args:
        tabular_data: [safe, shim, reg], where each is the increasing list of
            (%, $) returned by tabular_data().

    Returns:
        The filename of the final typeset PDF.
    """
    COLUMNS = 8
    safe, shim, reg = tabular_data
    # First half, rounded up to nearest multiple of COLUMNS
    regone = reg[:int(math.ceil(float(len(reg)) / (2 * COLUMNS))) * COLUMNS]
    # Everything from there on
    regtwo = reg[len(regone):]

    with open("template-worthtable.tex") as f:
        tex = f.read()

    tex = tex.replace("(safedata)", tex_one_table(safe, COLUMNS))
    tex = tex.replace("(shimdata)", tex_one_table(shim, COLUMNS))
    tex = tex.replace("(regdataone)", tex_one_table(regone, COLUMNS))
    tex = tex.replace("(regdatatwo)", tex_one_table(regtwo, COLUMNS))

    tex = tex.replace("(safeworth)", "\\${:.2f}".format(safe_worth))
    tex = tex.replace("(shimworth)", "\\${:.2f}".format(shim_worth))
    tex = tex.replace("(regworth)", "\\${:.2f}".format(reg_worth))  # twice

    filename = available_filename("worthtable")
    with open(filename, 'w') as f:
        f.write(tex)

    subprocess.call(["pdflatex", filename, "-interaction", "batchmode"])

    return filename.replace("tex", "pdf")

def tex_one_table(data, columns):
    """
    Generate the meat of the LaTeX tabular environment for one page of one rod
    (i.e. all of the safe or shim, or half of the reg).

    Args:
        data: [(%, $), ...] for one page of one rod.
        columns: The number of tuples wide the table should be.

    Returns:
        The rows of a LaTeX table 2*columns wide.
    """
    rows = int(math.ceil(float(len(data)) / columns))
    cols = [data[i:i + rows] for i in xrange(0, len(data), rows)]
    tex = ""
    for i in xrange(rows):
        tex += " & ".join("{0:.1f} & \\${1:.2f}".format(*col[i]) if
                          i < len(col) else " & " for col in cols)
        tex += "\\\\\n"  # that is, LaTeX's \\ and a hard newline
    return tex

def available_filename(prefix):
    """
    Return the first available filename of the form prefix.tex or prefix-N.tex
    where N is a positive integer.

    Args:
        prefix: How the filename should start.
    """
    if not os.path.exists("{}.tex".format(prefix)):
        return "{}.tex".format(prefix)

    i = 1
    while True:
        filename = "{}-{}.tex".format(prefix, i)
        if not os.path.exists(filename):
            return filename
        i += 1


def collect_rod(start):
    result = [start]
    while True:
        height = float(raw_input("Height after pull (percent): "))
        period = float(raw_input("                Period (ms): ")) / 1000.0
        result.append((height, period))
        print

        # Stop if we're full out (100% +/- tolerance defined by SOP)
        if height > 99.4:
            return result


def open_file(filename):
    if sys.platform.startswith("darwin"):
        subprocess.call(["open", filename])
    elif sys.platform.startswith("win"):
        os.startfile(filename)
    elif sys.platform.startswith("linux"):
        subprocess.call(["xdg-open", filename])
    # Couldn't identify the platform -- give up on opening the file.


def main():
    # Data entry stage
    print "This program will prompt you to enter data from SOP 34A (Control Rod"
    print "Calibration Form). You should already have filled it out, through"
    print "the bottom of Page 3. Not all the values on the form will be used."
    print
    print "Rod withdrawal times (seconds)"
    print "------------------------------"
    withdrawal_times = [float(raw_input(prompt)) for prompt in
                        "Safe: ", "Shim: ", " Reg: "]
    print

    print "Critical rod heights with Reg at bottom (percent)"
    print "-------------------------------------------------"
    critical_heights = [float(raw_input(prompt)) for prompt in
                        "Safe: ", "Shim: ", " Reg: "]
    print

    print "Reg rod"
    print "-------"
    reg = collect_rod(critical_heights[2])

    print "Safe and Shim rods"
    print "------------------"
    print "Calibrate togehter or separately?"
    print " [1] Together (this is the normal procedure)"
    print " [2] Separately"
    together = raw_input("Enter 1 or 2: ")
    print

    if together == "1":
        print "Critical rod heights with Shim and Reg at top (percent)"
        print "-------------------------------------------------------"
        critical_heights = [float(raw_input(prompt)) for prompt in
                            "Safe: ", "Shim: ", " Reg: "]
        print

        safe = [critical_heights[0]]
        shim = [critical_heights[1]]
        while True:
            safe_height = float(raw_input(" Safe height after pull (percent): "))
            period = (float(raw_input("                      Period (ms): ")) /
                      1000.0)
            shim_height = float(raw_input("Shim height when stable (percent): "))
            print

            safe.append((safe_height, period))
            shim[0] = (shim[0], period)
            shim.insert(0, shim_height)
            if safe_height > 99.4:
                break

    elif together == "2":
        print "Critical rod heights with Safe at bottom (percent)"
        print "--------------------------------------------------"
        critical_heights = [float(raw_input(prompt)) for prompt in
                            "Safe: ", "Shim: ", " Reg: "]
        print
        print "Safe rod"
        print "--------"
        safe = collect_rod(critical_heights[0])

        print "Critical rod heights with Shim at bottom (percent)"
        print "--------------------------------------------------"
        critical_heights = [float(raw_input(prompt)) for prompt in
                            "Safe: ", "Shim: ", " Reg: "]
        print
        print "Shim rod"
        print "--------"
        shim = collect_rod(critical_heights[1])

    else:
        return

    rod_pulls = [safe, shim, reg]

    # Done with data entry: start calculating.
    reactivity_results = [reactivity_calc(rod[0], rod[1:]) for rod in rod_pulls]
    [(safe_int, safe_int_data, safe_add, safe_add_data),
     (shim_int, shim_int_data, shim_add, shim_add_data),
     (reg_int, reg_int_data, reg_add, reg_add_data)] = reactivity_results

    for rod, result in zip(["safe", "shim", "reg"], reactivity_results):
        plot(*((rod,) + result))

    safe_full = safe_int(0)
    shim_full = shim_int(0)
    reg_full = reg_int(0)
    table_filename = tex_tables(tabular_data([safe_int, shim_int, reg_int]),
                                safe_full, shim_full, reg_full)
    print "Filename: " + table_filename
    open_file(table_filename)

    report_filename = tex_report(tech_specs(safe_int, safe_add, shim_int,
                                            shim_add, reg_int, reg_add,
                                            rod_pulls[0][0], rod_pulls[1][0],
                                            withdrawal_times))
    print "Filename: " + report_filename
    open_file(report_filename)

if __name__ == "__main__":
    main()
