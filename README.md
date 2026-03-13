# Numerical Root Finding Toolbox

## Project Description

This project implements several classical numerical methods used to solve nonlinear equations of the form:

f(x) = 0

The toolbox is implemented in Python using only standard Python libraries and the `cmath` module for complex numbers.

All algorithms are implemented manually without using external numerical libraries such as NumPy or SciPy.

The project provides a reusable class called `RootFindingProblem` which solves equations using different numerical algorithms through one public interface called `solve()`.

---

## Implemented Methods

The following methods are implemented:

- Bisection Method
- Fixed Point Iteration
- Newton’s Method
- Secant Method
- False Position Method (Regula Falsi)
- Horner’s Method
- Muller’s Method

---

## Algorithms Overview

### Bisection Method
The bisection method repeatedly divides an interval in half and selects the subinterval where the function changes sign.

### Fixed Point Iteration
This method rewrites the equation in the form:

x = g(x)

and repeatedly computes the next approximation until convergence.

### Newton’s Method
Newton’s method uses the derivative of the function:

x(n+1) = x(n) - f(x)/f'(x)

to approximate the root.

### Secant Method
The secant method approximates the derivative using two previous points instead of computing the analytical derivative.

### False Position Method
This method improves the bisection method by using linear interpolation between two points.

### Horner’s Method
Horner’s method efficiently evaluates polynomials with fewer multiplications.

### Muller’s Method
Muller’s method uses quadratic interpolation with three points and can find both real and complex roots.

---

## Project Structure
```
root-finding-project/

root_find.py
examples.py
README.md
```

## How to Run the Examples
```bash
python examples.py
```
---

## Example Usage

Example of using Newton's method:

```python
from root_find import RootFindingProblem

def f(x):
    return x**3 - x - 2

def df(x):
    return 3*x**2 - 1

p = RootFindingProblem(f=f, df=df)

root = p.solve("newton", x0=1.5)

print("Root:", root)
```
---

## Requirements

- Python 3
- Standard Python libraries
- cmath module for complex numbers

No external numerical libraries are used.

---

## Author

This project was created for a Numerical Analysis mini-project assignment.
