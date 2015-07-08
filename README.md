# CD\*

CD\* is an expression compiler written in C++ that generates C code of functions that contain derivatives computation.  For example, the following code:
```
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;
    auto xCu = xSq * x;
    // y = 2x^3 + 3x^2 + 4x + 5
    auto y = 2.0 * xCu + 3.0 * xSq + 4.0 * x + 5.0; 
    auto dydx = std::make_shared<NamedAssignment>("dydx", 0, Derivatives({{y, x}})[0]);
    EmitFunction({xArg}, {dydx}, "derv", std::cout);
```
Generates the following C code in standard output
```
void derv(const double x, double dydx[1]) {
    double t[7];
    t[0] = x * x;
    t[1] = 2.000000 * t[0];
    t[2] = 2.000000 * x;
    t[3] = 3.000000 + t[2];
    t[4] = t[3] * t[2];
    t[5] = t[1] + t[4];
    t[6] = 4.000000 + t[5];
    dydx[0] = t[6];
}
```
It is possible to generate higher-order derivatives or mixed derivatives by mixing the calls of the Derivatives function.  See test.cpp for more examples.

CD\* implements the [D\* algorithm](http://dl.acm.org/citation.cfm?id=1276512) to efficiently compute the derivatives of a function.  The algorithm is efficient for small-to-medium dimensional (say, <100) functions because it performs costly path-finding algorithm on the expression graph.  If you need to generate the code of the Jacobian or Hessian of a function with 10000  variables, you should also take a look at the following libraries: [CasADi](https://github.com/casadi/casadi/wiki), [theano](http://deeplearning.net/software/theano/), [CppADCode](https://github.com/joaoleal/CppADCodeGen/)

This project is still in its very early stage, the interface may change significantly in future development.

More to come.
