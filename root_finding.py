import cmath


class RootFindingProblem:

    def __init__(self, f=None, df=None, g=None):
        self.f = f
        self.df = df
        self.g = g


    def solve(self, method, **kwargs):

        if method == "bisection":
            return self._bisection(**kwargs)

        elif method == "fixed_point":
            return self._fixed_point(**kwargs)

        elif method == "newton":
            return self._newton(**kwargs)

        elif method == "secant":
            return self._secant(**kwargs)

        elif method == "false_position":
            return self._false_position(**kwargs)

        elif method == "muller":
            return self._muller(**kwargs)

        else:
            raise ValueError("Unknown solving method")


# --------------------------------------------------
# Bisection Method
# --------------------------------------------------

    def _bisection(self, a, b, tol=1e-6, max_iter=100):

        if self.f(a) * self.f(b) > 0:
            raise ValueError("Invalid interval for bisection")

        for _ in range(max_iter):

            c = (a + b) / 2

            if abs(self.f(c)) < tol or abs(b - a) < tol:
                return c

            if self.f(a) * self.f(c) < 0:
                b = c
            else:
                a = c

        raise ValueError("No convergence after max_iter")


# --------------------------------------------------
# Fixed Point Iteration
# --------------------------------------------------

    def _fixed_point(self, x0, tol=1e-6, max_iter=100):

        if self.g is None:
            raise ValueError("Missing fixed-point function g(x)")

        x = x0

        for _ in range(max_iter):

            x_new = self.g(x)

            if abs(x_new - x) < tol:
                return x_new

            x = x_new

        raise ValueError("No convergence after max_iter")


# --------------------------------------------------
# Newton Method
# --------------------------------------------------

    def _newton(self, x0, tol=1e-6, max_iter=100):

        if self.df is None:
            raise ValueError("Missing derivative df(x) for Newton method")

        x = x0

        for _ in range(max_iter):

            d = self.df(x)

            if d == 0:
                raise ZeroDivisionError("Division by zero in Newton method")

            x_new = x - self.f(x) / d

            if abs(x_new - x) < tol:
                return x_new

            x = x_new

        raise ValueError("No convergence after max_iter")


# --------------------------------------------------
# Secant Method
# --------------------------------------------------

    def _secant(self, x0, x1, tol=1e-6, max_iter=100):

        for _ in range(max_iter):

            f0 = self.f(x0)
            f1 = self.f(x1)

            if (f1 - f0) == 0:
                raise ZeroDivisionError("Division by zero in Secant method")

            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

            if abs(x2 - x1) < tol:
                return x2

            x0 = x1
            x1 = x2

        raise ValueError("No convergence after max_iter")


# --------------------------------------------------
# False Position Method
# --------------------------------------------------

    def _false_position(self, a, b, tol=1e-6, max_iter=100):

        if self.f(a) * self.f(b) > 0:
            raise ValueError("Invalid interval for False Position")

        for _ in range(max_iter):

            fa = self.f(a)
            fb = self.f(b)

            if (fb - fa) == 0:
                raise ZeroDivisionError("Division by zero in False Position")

            c = (a * fb - b * fa) / (fb - fa)

            if abs(self.f(c)) < tol:
                return c

            if fa * self.f(c) < 0:
                b = c
            else:
                a = c

        raise ValueError("No convergence after max_iter")


# --------------------------------------------------
# Horner Method (Polynomial Evaluation)
# --------------------------------------------------

    def _horner(self, coeffs, x):

        result = coeffs[0]

        for c in coeffs[1:]:
            result = result * x + c

        return result


# --------------------------------------------------
# Muller Method (Complex roots allowed)
# --------------------------------------------------

    def _muller(self, x0, x1, x2, tol=1e-6, max_iter=100):

        for _ in range(max_iter):

            f0 = self.f(x0)
            f1 = self.f(x1)
            f2 = self.f(x2)

            h0 = x1 - x0
            h1 = x2 - x1

            if h0 == 0 or h1 == 0:
                raise ZeroDivisionError("Invalid starting points for Muller")

            d0 = (f1 - f0) / h0
            d1 = (f2 - f1) / h1

            if (h1 + h0) == 0:
                raise ZeroDivisionError("Muller denominator became zero")

            a = (d1 - d0) / (h1 + h0)
            b = a * h1 + d1
            c = f2

            disc = cmath.sqrt(b*b - 4*a*c)

            if abs(b + disc) > abs(b - disc):
                denom = b + disc
            else:
                denom = b - disc

            if denom == 0:
                raise ZeroDivisionError("Division by zero in Muller method")

            dx = -2 * c / denom

            x3 = x2 + dx

            if abs(dx) < tol:
                return x3

            x0 = x1
            x1 = x2
            x2 = x3

        raise ValueError("No convergence after max_iter")