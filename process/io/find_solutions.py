"""Tests for the Process solver."""
import time
from numpy import histogram
import process
from process.io.mfile import MFile
from process.io.python_fortran_dicts import get_dicts
from process.io.process_config import TestProcessConfig
from process.io.process_funcs import (
    get_neqns_itervars,
    update_ixc_bounds,
    get_variable_range,
    check_input_error,
    vary_iteration_variables,
    process_stopped,
    get_solution_from_mfile,
    process_warnings,
)
from pathlib import Path

# Load dicts from dicts JSON file
process_dicts = get_dicts()
IFAIL_SUCCESS = process_dicts["IFAIL_SUCCESS"]


def find_solutions(config_path, solver_config="solver.toml"):
    """Code to test PROCESS Solver by choosing different starting values
    for the iteration parameters around the initial values in the INPUT file
    Code modified by Sarah Medley in April 2015 to also calculate Q (fusion power/
    injected power) and output this to SolverTest.out

    Author: H. Lux (Hanni.Lux@ccfe.ac.uk)

    Date: November 2013 (Initial version)
        March 2014 (Initial release version)
        April 2015 - Code modified by Sarah Medley to also calculate Q - see above

    - Input files -
    test_process.conf: in integration/data dir

    - Output files -
    Saved to temporary test dir
    OUT.DAT     - PROCESS output of last run
    PLOT.DAT    - PROCESS output of last run
    MFILE.DAT   - PROCESS output of last run
    process.log - logfile of PROCESS output to stdout
    time.info   - stores the run time of the inner loop
    SolverTest.out -

    Ifail values:
    -1: Error in logfile, no solution vector in OUT.DAT
    0: Error in VMCON: improper input parameters
    1: normal run
    2: too many function calls in VMCON
    3/4: either constraint/objective function + gradients
        are not calculated well or data too noisy
    5: either no feasible solution to the problem,
        or identity matrix is a bad approximation ot the Hessian
    6: restriction by an artificial bound or singular matrix
        in quadratic problem.

    :param wdir: temp test dir containing integration test data
    :type wdir: Path
    """
    converged_runs_count = 0

    # Path to varied input file
    input_path = Path("IN.DAT")
    CONFIG = TestProcessConfig(config_path)
    CONFIG.setup()

    process.fortran.init_module.init_all_module_vars()
    process.fortran.init_module.init()

    NEQNS, ITERVARS = get_neqns_itervars()

    update_ixc_bounds()

    LBS, UBS = get_variable_range(ITERVARS, CONFIG.factor)

    TABLE_ERR = []
    TABLE_OBJ = []
    TABLE_CON = []
    TABLE_IN = []
    TABLE_OUT = []
    TABLE_POWF = []
    TABLE_PINJ = []
    # TABLE_RATIO_POWF_PINJ = []

    WARNING_CNT = 0

    START_TIME = time.time()

    for i in range(CONFIG.niter):
        # modify IN.DAT each time
        input_values = vary_iteration_variables(ITERVARS, LBS, UBS)

        TABLE_IN += [input_values]

        print(i, end=" ")

        # Actually perform the Process run for this varied input file
        CONFIG.run_process(input_path)

        check_input_error()

        if process_stopped():
            # Probably a model exception
            ifail, objective_function, constraints, table_sol, table_res = (
                -1,
                "0",
                "0",
                ["0"] * len(ITERVARS),
                ["0"] * NEQNS,
            )
            fusion_power = "0"
            injected_power = "0"
            # ratio_fus_inj_power = 0
        else:
            # Solver returned; could have converged or not
            (
                ifail,
                objective_function,
                constraints,
                table_sol,
                table_res,
            ) = get_solution_from_mfile(NEQNS, len(ITERVARS))

            # Access MFILE.DAT
            m_file = MFile(filename="MFILE.DAT")
            ind = -1  # last scan point

            # Obtain fusion power and injected power from MFILE.DAT
            fusion_power = m_file.data["powfmw"].get_scan(ind)
            injected_power = m_file.data["pinjmw"].get_scan(ind)
            # Removed: can cause div by 0 errors
            # ratio_fus_inj_power = fusion_power / injected_power

            if ifail == IFAIL_SUCCESS and process_warnings():
                WARNING_CNT += 1

            if ifail == 1:
                # Process converged: save the mfile
                print(f"Process converged. {converged_runs_count = }")
                current_mfile = Path("MFILE.DAT")
                renamed_mfile = Path(f"sol_{converged_runs_count}_MFILE.DAT")
                current_mfile.rename(renamed_mfile)
                converged_runs_count += 1

        TABLE_OUT += [table_sol + table_res]
        TABLE_ERR += [ifail]
        TABLE_OBJ += [objective_function]
        TABLE_CON += [constraints]
        TABLE_POWF += [fusion_power]
        TABLE_PINJ += [injected_power]
        # TABLE_RATIO_POWF_PINJ += [ratio_fus_inj_power]

        i += 1

    END_TIME = time.time()

    ##########
    TIMEFILE = open("time.info", "w")
    TIMEFILE.write("wall clock time of run %f s" % (END_TIME - START_TIME))
    TIMEFILE.close()

    #########
    print("\n Error codes:")
    ERRHIST, ERRBINS = histogram(
        TABLE_ERR, bins=[-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    )

    for i in range(len(ERRHIST)):
        print(
            "code %2i: %4i out of %4i = %3.1f "
            % (
                int((ERRBINS[i + 1] + ERRBINS[i]) / 2.0),
                ERRHIST[i],
                CONFIG.niter,
                ERRHIST[i] * 100.0 / CONFIG.niter,
            )
            + "%"
        )

    if WARNING_CNT > 0:
        # TODO: calculate fraction of warnings in sucessful runs.
        print("\nIn %f%% of the successful runs warnings occurred." % WARNING_CNT)

    ##########
    OUTFILE = open("SolverTest.out", "w")

    OUTHEADER = "#Err\t ObjectiveF\t Constraints\t"
    for inp_str in ITERVARS:  # input variables
        OUTHEADER += "%s_in\t" % inp_str
    for var_no in range(len(ITERVARS)):  # output variables
        OUTHEADER += "itvar%03i\t" % (var_no + 1)
    for con_no in range(NEQNS):  # constraints
        OUTHEADER += "constr%03i\t" % (con_no + 1)
    # OUTHEADER += 'Fusion_Power_(MW)\t'
    # OUTHEADER += 'Injected_Power_(MW)\t'
    # OUTHEADER += "Ratio_fus_inj_power_(MW)"
    OUTHEADER += "\n"
    OUTFILE.write(OUTHEADER)

    for i in range(CONFIG.niter):
        outstr = "%d\t %s\t %s\t" % (TABLE_ERR[i], TABLE_OBJ[i], TABLE_CON[i])
        for valuein in TABLE_IN[i]:
            outstr += "%s\t" % valuein
        for valueout in TABLE_OUT[i]:
            outstr += "%s\t" % valueout
        # outstr += '%d\t' % TABLE_POWF[i]
        # outstr += '%d\t' % TABLE_PINJ[i]
        # outstr += "%.1f\t" % TABLE_RATIO_POWF_PINJ[i]
        outstr += "\n"
        OUTFILE.write(outstr)

    OUTFILE.close()
