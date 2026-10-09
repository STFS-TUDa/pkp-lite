"""Generic ODE reactor: integrates dY/dt = rate_fn(T, Y) along a piecewise-
linear T(t) profile (the `operating_point`)."""

import numpy as np
from scipy.integrate import solve_ivp

# Fitting explores unphysical coefficients that can produce degenerate ODE
# steps (e.g. 0/0 in DOP853); harmless, so silence the console spam.
np.seterr(divide='ignore', invalid='ignore')

class reactor:

    def __init__(self,fuel,rate_function,rate_function_arguments=None):

        self.fuel = fuel
        self.rate_function = rate_function
        self.rate_function_arguments=rate_function_arguments

    def run(self,operating_point):
        """Integrate the state Y (starting at `fuel`) over one operating point.

        operating_point: list of [t, T] points defining a piecewise-linear
        temperature profile; T is appended to Y as an extra state so the
        ODE solver advances both together.
        """

        self.operating_point = np.array(operating_point)

        self.dTdt_array = (np.diff(self.operating_point[:, 1]) /
                            np.diff(self.operating_point[:, 0]))

        # is used in the calc_rate function
        self.index = 0

        y0 = np.hstack((self.fuel, self.operating_point[0,1]))
        t0 = self.operating_point[0,0]
        tend = self.operating_point[-1,0]
        t_eval = np.arange(t0,tend+1e-4,np.minimum((tend-t0)/100, 1e-4))

        # Wrapper for solve_ivp
        def rhs(t, y):
            return self.calc_rate(t, y)

        # Integrate
        sol = solve_ivp(rhs, (t0, tend), y0, method="DOP853", dense_output=True, t_eval=t_eval, max_step=1e-4)

        # Reconstruct result array in the same style you had
        # include both times and state values
        result = np.column_stack((sol.t, sol.y.T))

        return result

    def calc_rate(self,t,y):

        # decompose T and y 
        T = y[-1]
        Y = y[:-1]

        # calculate rate
        rate = self.rate_function(T, Y,**self.rate_function_arguments)

        # advance dTdt if necessary
        if T > self.operating_point[self.index + 1,1]:
            self.index = int(np.minimum(self.index+1, len(self.dTdt_array)-1))

        dTdt = self.dTdt_array[self.index]
        # construct rate
        y_rate = np.hstack((rate,dTdt))

        return y_rate