from root_finding import RootFindingProblem


# -----------------------------------------
# Example 1 : Real function for most methods
# f(x) = x^3 - x - 2
# -----------------------------------------

def f(x):
    return x**3 - x - 2


def df(x):
    return 3*x**2 - 1


def g(x):
    return (x + 2)**(1/3)


p = RootFindingProblem(f=f, df=df, g=g)


print("Bisection:", p.solve("bisection", a=1, b=2))

print("Newton:", p.solve("newton", x0=1.5))

print("Secant:", p.solve("secant", x0=1, x1=2))

print("False Position:", p.solve("false_position", a=1, b=2))

print("Fixed Point:", p.solve("fixed_point", x0=1))


# -----------------------------------------
# Example 2 : Muller Method (complex root)
# f(x) = x^2 + 1
# roots = i and -i
# -----------------------------------------

def f_complex(x):
    return x**2 + 1


m = RootFindingProblem(f=f_complex)

print("Muller:", m.solve("muller", x0=0, x1=1, x2=2))


# -----------------------------------------
# Example 3 : Horner Method
# Polynomial: x^3 - 6x^2 + 11x - 6
# -----------------------------------------

coeffs = [1, -6, 11, -6]

h = RootFindingProblem()

print("Horner evaluation at x=2:", h._horner(coeffs, 2))