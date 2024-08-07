"""An adapter for different solvers."""

import logging
from process.fortran import numerics, global_variables
from process.utilities.f2py_string_patch import f2py_compatible_to_string
import numpy as np
from process.evaluators import Evaluators
from abc import ABC, abstractmethod
from typing import Optional, Union
import importlib
from pyvmcon import (
    AbstractProblem,
    Result,
    solve,
    QSPSolverException,
    VMCONConvergenceException,
    LineSearchConvergenceException,
)
import nlopt
import time

logger = logging.getLogger(__name__)


class _Solver(ABC):
    """Base class for different solver implementations.

    :param ABC: abstract base class
    :type ABC: ABC
    """

    def __init__(self) -> None:
        """Initialise a solver."""
        # Exit code for the solver
        self.ifail = 0
        self.tolerance = numerics.epsvmc
        self.b: Union[float, None] = None

    def set_evaluators(self, evaluators: Evaluators) -> None:
        """Set objective and constraint functions and their gradient evaluators.

        :param evaluators: objective and constraint evaluators
        :type evaluators: Evaluators
        """
        self.evaluators = evaluators

    def set_opt_params(self, x_0: np.ndarray) -> None:
        """Define the initial optimisation parameters.

        :param x_0: optimisation parameters vector
        :type x_0: np.ndarray
        """
        self.x_0 = x_0

    def set_bounds(
        self,
        bndl: np.ndarray,
        bndu: np.ndarray,
        ilower: Optional[np.ndarray] = None,
        iupper: Optional[np.ndarray] = None,
    ) -> None:
        """Set the bounds on the optimisation parameters.

        :param bndl: lower bounds for the optimisation parameters
        :type bndl: np.ndarray
        :param bndu: upper bounds for the optimisation parameters
        :type bndu: np.ndarray
        :param ilower: array of 0s and 1s to activate lower bounds on
        optimsation parameters in x
        :type ilower: np.ndarray, optional
        :param iupper: array of 0s and 1s to activate upper bounds on
        optimsation parameters in x
        :type iupper: np.ndarray, optional
        """
        self.bndl = bndl
        self.bndu = bndu

        # TODO Remove ilower/iupper and use finite vs. infinite values in bndl/bndu
        # instead to determine defined vs. undefined bounds
        # If lower and upper bounds switch arrays aren't specified, set all bounds
        # to defined
        if ilower is None and iupper is None:
            ilower = np.ones(bndl.shape[0], dtype=int)
            iupper = np.ones(bndu.shape[0], dtype=int)

        self.ilower = ilower
        self.iupper = iupper

    def set_constraints(self, m: int, meq: int) -> None:
        """Set the total number of constraints and equality constraints.

        :param m: number of constraint equations
        :type m: int
        :param meq: of the constraint equations, how many are equalities
        :type meq: int
        """
        self.m = m
        self.meq = meq

    def set_tolerance(self, tolerance: float) -> None:
        """Set tolerance for solver termination.

        :param tolerance: tolerance for solver termination
        :type tolerance: float
        """
        self.tolerance = tolerance

    def set_b(self, b: float) -> None:
        """Set the multiplier for the Hessian approximation.

        :param b: multiplier for an identity matrix as input for the Hessian b(n,n)
        :type b: float
        """
        self.b = b

    @abstractmethod
    def solve(self) -> int:
        """Run the optimisation.

        :return: solver error code
        :rtype: int
        """
        pass


class VmconProblem(AbstractProblem):
    def __init__(self, evaluator, nequality, ninequality) -> None:
        self._evaluator = evaluator
        self._nequality = nequality
        self._ninequality = ninequality

    def __call__(self, x: np.ndarray) -> Result:
        n = x.shape[0]
        objf, conf = self._evaluator.fcnvmc1(n, self.total_constraints, x, 0)
        fgrd, cnorm = self._evaluator.fcnvmc2(n, self.total_constraints, x, n)

        return Result(
            objf,
            fgrd,
            conf[: self.num_equality],
            cnorm[:, : self.num_equality].T,
            conf[self.num_equality :],
            cnorm[:, self.num_equality :].T,
        )

    @property
    def num_equality(self) -> int:
        return self._nequality

    @property
    def num_inequality(self) -> int:
        return self._ninequality


class Vmcon(_Solver):
    """New VMCON implementation.

    :param _Solver: Solver base class
    :type _Solver: _Solver
    """

    def solve(self) -> int:
        """Optimise using new VMCON.

        :return: solver error code
        :rtype: int
        """
        problem = VmconProblem(self.evaluators, self.meq, self.m - self.meq)

        B = None
        if self.b is not None:
            B = np.identity(numerics.nvar) * self.b

        def _solver_callback(i: int, _result, _x, convergence_param: float):
            numerics.nviter = i + 1
            global_variables.convergence_parameter = convergence_param
            print(
                f"{i+1} | Convergence Parameter: {convergence_param:.3E}",
                end="\r",
                flush=True,
            )

        def _ineq_cons_satisfied(
            result: Result,
            _x: np.ndarray,
            _delta: np.ndarray,
            _lambda_eq: np.ndarray,
            _lambda_in: np.ndarray,
        ) -> bool:
            """Check that inequality constraints are satisfied.

            This additional convergence criterion ensures that solutions have
            satisfied inequality constraints.

            :param result: evaluation of current optimisation parameter vector
            :type result: Result
            :param _x: current optimisation parameter vector
            :type _x: np.ndarray
            :param _delta: search direction for line search
            :type _delta: np.ndarray
            :param _lambda_eq: equality Lagrange multipliers
            :type _lambda_eq: np.ndarray
            :param _lambda_in: inequality Lagrange multipliers
            :type _lambda_in: np.ndarray
            :return: True if inequality constraints satisfied
            :rtype: bool
            """
            # Check all ineqs positive, i.e. satisfied
            if np.all(result.ie >= 0.0):
                return True
            else:
                return False

        try:
            x, _, _, res = solve(
                problem,
                self.x_0,
                self.bndl,
                self.bndu,
                max_iter=global_variables.maxcal,
                epsilon=self.tolerance,
                qsp_options={"eps_rel": 1e-1, "adaptive_rho_interval": 25},
                initial_B=B,
                callback=_solver_callback,
                additional_convergence=_ineq_cons_satisfied,
            )
        except VMCONConvergenceException as e:
            if isinstance(e, LineSearchConvergenceException):
                self.info = 3
            elif isinstance(e, QSPSolverException):
                self.info = 5
            else:
                self.info = 2

            logger.warning(str(e))

            x = e.x
            res = e.result

        except ValueError as e:
            itervar_name_list = ""
            for count, iter_var in enumerate(numerics.ixc[: numerics.nvar]):
                itervar_name = f2py_compatible_to_string(numerics.lablxc[iter_var - 1])
                itervar_name_list += f"{count}: {itervar_name} \n"

            logger.warning(f"Active iteration variables are : \n{itervar_name_list}")
            raise e

        else:
            self.info = 1

        # print a blank line because of the carridge return
        # in the callback
        print()

        self.x = x
        self.objf = res.f
        self.conf = np.hstack((res.eq, res.ie))

        return self.info


class VmconBounded(Vmcon):
    """A solver that uses VMCON but checks x is in bounds before running"""

    def set_opt_params(self, x_0: np.ndarray) -> None:
        lower_violated = np.less(x_0, self.bndl)
        upper_violated = np.greater(x_0, self.bndu)

        for index, entry in enumerate(lower_violated):
            if entry:
                x_0[index] = self.bndl[index]
        for index, entry in enumerate(upper_violated):
            if entry:
                x_0[index] = self.bndu[index]
        self.x_0 = x_0


class SLSQP(_Solver):
    def obj_func(self, x, grad):
        # Must be passed these args from nlopt requirements
        objf, conf = self.evaluators.fcnvmc1(self.n, self.m, x, self.ifail)
        self.constr_res = np.sqrt(
            np.sum(np.square(conf[: self.meq]))
        )  # only for equality constraints

        if grad.size > 0:
            # Gradient required by solver; modify grad in-place
            fgrd, cnorm = self.evaluators.fcnvmc2(self.n, self.m, x, self.lcnorm)
            grad[...] = fgrd

        self.eval_count += 1
        status = (
            f"Evaluation {self.eval_count}, objective function = {objf:.5}, "
            f"constraint residuals = {self.constr_res:.3e}"
        )
        print(status)
        logger.info(status)
        return objf

    def constraint_eq_vec(self, result, x, grad):
        # i is no. of constraint
        objf, conf = self.evaluators.fcnvmc1(self.n, self.m, x, self.ifail)
        self.constr_res = np.sqrt(
            np.sum(np.square(conf[: self.meq]))
        )  # only for equality constraints

        # Check if constraints are below tolerance yet
        # conf_gt_tol = np.sum((np.abs(conf[:meq]) > CONSTR_TOL))
        # if conf_gt_tol == 0:
        #     logger.info("Constraints are all below tolerance")
        # else:
        #     logger.info(f"Constraints above tolerance: {conf_gt_tol=}")

        if grad.size > 0:
            # Gradient required by solver; modify grad in-place
            fgrd, cnorm = self.evaluators.fcnvmc2(self.n, self.m, x, self.lcnorm)
            if cnorm.ndim == 1:
                # 1 constraint, 1 optimisation parameter (used in tests)
                # TODO Add np.newaxis to fcnvmc2 to avoid this
                grad[...] = -cnorm[:]
            else:
                grad[...] = np.transpose(-cnorm[: x.shape[0], : self.meq])

        # Negative conf and cnorm: nlopt requires opposite formulation of
        # VMCON's (and Process's) constraint form
        result[...] = -conf[: self.meq]

    def constraint_ineq_vec(self, result, x, grad):
        objf, conf = self.evaluators.fcnvmc1(self.n, self.m, x, self.ifail)
        if grad.size > 0:
            # Gradient required by solver; modify grad in-place
            fgrd, cnorm = self.evaluators.fcnvmc2(self.n, self.m, x, self.lcnorm)
            if cnorm.ndim == 1:
                # 1 constraint, 1 optimisation parameter (used in tests)
                # TODO Add np.newaxis to fcnvmc2 to avoid this
                grad[...] = -cnorm[0]
            else:
                grad[...] = np.transpose(-cnorm[: x.shape[0], self.meq : self.m])

        # Negative conf and cnorm: nlopt requires opposite formulation of
        # VMCON's (and Process's) constraint form
        result[...] = -conf[self.meq : self.m]

    def solve(self) -> int:
        """Try running nlopt."""
        # global objf, conf, fgrd, cnorm
        self.eval_count = 0
        self.constr_res = 0.0

        self.n = self.x_0.shape[0]
        self.lcnorm = 176  # (ippnvars + 1), but could just be n!

        # Solver tolerances (high)
        MAIN_TOL = 1e-8
        # LOCAL_TOL = 1e-8 # for AUGLAG method only
        CONSTR_TOL = 1e-10

        # Set up optimiser object
        # opt = nlopt.opt(nlopt.AUGLAG, n)
        opt = nlopt.opt(nlopt.LD_SLSQP, self.n)

        # Set subsidiary optimiser
        # For AUGLAG only
        # local_opt = nlopt.opt(nlopt.LD_SLSQP, n)

        opt_name = opt.get_algorithm_name()
        # local_opt_name = local_opt.get_algorithm_name()
        # logger.info(f"{opt_name=}, {local_opt_name=}")
        logger.info(f"{opt_name=}")
        logger.info(f"{MAIN_TOL=:.3e}, {CONSTR_TOL=:.3e}")
        # logger.info(f"{MAIN_TOL=:.3e}, {LOCAL_TOL=:.3e}, {CONSTR_TOL=:.3e}")

        # Define tolerances
        opt.set_ftol_rel(MAIN_TOL)
        # local_opt.set_ftol_rel(LOCAL_TOL)

        # Need to terminate in solver test case 3!
        opt.set_maxeval(1000)

        # opt.set_local_optimizer(local_opt)

        # if self.maximise:
        #     # Maximisation required for test case 5
        #     opt.set_max_objective(self.obj_func)
        # else:
        opt.set_min_objective(self.obj_func)

        # Check bounds are activated for all optimisation parameters (default case)
        # If not, handle it
        if not (np.all(self.ilower) and np.all(self.iupper)):
            # Some bounds are inactive: used e.g. in solver integration tests
            for i in range(self.ilower.shape[0]):
                if self.ilower[i] == 0:
                    # Inactive lower bound: set to -inf
                    self.bndl[i] = -np.inf

            for i in range(self.iupper.shape[0]):
                if self.iupper[i] == 0:
                    # Inactive upper bound: set to +inf
                    self.bndu[i] = np.inf

        opt.set_lower_bounds(self.bndl)
        opt.set_upper_bounds(self.bndu)

        # Kludge initial normalised x into the normalised bounds range if required
        for i in range(self.x_0.shape[0]):
            if self.x_0[i] < self.bndl[i]:
                self.x_0[i] = self.bndl[i]
            elif self.x_0[i] > self.bndu[i]:
                self.x_0[i] = self.bndu[i]

        # Constraints
        if self.meq > 0:
            eq_constr_tols = np.full(self.meq, CONSTR_TOL)
            opt.add_equality_mconstraint(self.constraint_eq_vec, eq_constr_tols)
        if self.meq < self.m:
            ineq_constr_tols = np.full((self.m - self.meq), CONSTR_TOL)
            opt.add_inequality_mconstraint(self.constraint_ineq_vec, ineq_constr_tols)

        start_time = time.time()
        x_opt = opt.optimize(self.x_0)
        end_time = time.time()
        duration = end_time - start_time

        opt_val = opt.last_optimum_value()
        return_value = opt.last_optimize_result()

        # Main opt iterations
        main_opt_evals = opt.get_numevals()

        if return_value < 0:
            raise RuntimeError("nlopt didn't converge")
        elif return_value == nlopt.MAXEVAL_REACHED:
            # For giving up in int test 3
            # TODO Reconcile MAXEVAL with above exception: int case 3 needs to pass
            info = 5
        elif return_value > 0:
            # TODO Might want to be aware of other return values
            info = 1

        print(f"{return_value=}")

        # Recalculate conf at optimum x
        objf, conf = self.evaluators.fcnvmc1(self.n, self.m, x_opt, self.ifail)
        self.constr_res = np.sqrt(
            np.sum(np.square(conf[: self.meq]))
        )  # only for equality constraints

        logger.info(f"Main opt evaluations = {main_opt_evals}")
        logger.info(f"{opt_val=:.3e}, {self.constr_res=:.3e}")
        logger.info(f"{duration=:.1f}\n")
        logger.info(f"{conf[:self.meq]=}")

        print(f"{self.constr_res=:.3e}")

        # Check how many conf elements are above tolerance
        conf_gt_tol = np.sum((np.abs(conf[: self.meq]) > CONSTR_TOL))
        logger.info(f"Constraints above tolerance: {conf_gt_tol}")

        # Store required results on object
        self.objf = objf
        self.conf = conf
        self.x = x_opt

        return info


def get_solver(solver_name: str = "vmcon") -> _Solver:
    """Return a solver instance.

    :param solver_name: solver to create, defaults to "vmcon"
    :type solver_name: str, optional
    :return: solver to use for optimisation
    :rtype: _Solver
    """
    solver: _Solver

    if solver_name == "vmcon":
        solver = Vmcon()
    elif solver_name == "vmcon_bounded":
        solver = VmconBounded()
    elif solver_name == "slsqp":
        solver = SLSQP()
    else:
        try:
            solver = load_external_solver(solver_name)
        except Exception as e:
            raise ValueError(
                f'Solver name is not an inbuilt PROCESS solver or recognised package "{solver_name}"'
            ) from e

    return solver


def load_external_solver(package: str):
    """Attempts to load a package of name `package`.

    If a package of the name is available, return the `__process_solver__`
    attribute of that package or raise an `AttributeError`."""
    module = importlib.import_module(package)

    solver = getattr(module, "__process_solver__", None)

    if solver is None:
        raise AttributeError(
            f"Module {module.__name__} does not have a '__process_solver__' attribute."
        )

    return solver()
