#include "cdstar/expression.h"
#include "cdstar/library.h"

#include <iostream>
#include <cmath>

using namespace cdstar;

void TestConstant() {
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto y = std::make_shared<NamedAssignment>("y", 0, std::make_shared<Constant>(1.0));
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_constant");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";
    lib.AddFunction(input, dervOutput, dervFuncName);         
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {1.0};
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestConstant Failed");
        }
    }   
    {
        double x = 0.5, y[1];
        df(x, y);
        double ref[1] = {0.0};        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestConstant Failed");
        }                        
    }       
}

void TestLinear() {
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto y = std::make_shared<NamedAssignment>("y", 0, x);
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_linear");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";
    lib.AddFunction(input, dervOutput, dervFuncName);         
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {x};
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestLinear Failed");
        }
    }
    {
        double x = 0.5, y[1];
        df(x, y);
        double ref[1] = {1.0};        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestLinear Failed");
        }       
    }   
}

void TestPolynomial() {
    // y = 2x^3 - 3x^2 + 4x + 5
    // dy/dx = 6x^2 - 6x + 4
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;
    auto xCu = xSq * x;
    auto y = std::make_shared<NamedAssignment>("y", 0, 2.0 * xCu - 3.0 * xSq + 4.0 * x + 5.0);    
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_polynomial");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {
            2.0 * x * x * x - 3.0 * x * x + 4.0 * x + 5.0
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestPolynomial Failed");
        }
    }   
    {
        double x = 0.5, y[1];
        df(x, y);
        double ref[1] = {
            6.0 * x * x - 6.0 * x + 4.0
        };
        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestPolynomial Failed");
        }                        
    }       
}

void TestRational() {
    // y = (x^2 + 2x + 3) / (4x^2 + 5x + 6)
    // dy/dx = - 3 * (x^2 + 4x + 1) / (4x^2 + 5x + 6)^2
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;
    auto y = std::make_shared<NamedAssignment>("y", 0, (xSq + 2.0 * x + 3.0) / (4.0 * xSq + 5.0 * x + 6.0));
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_rational");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";
    lib.AddFunction(input, dervOutput, dervFuncName);
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {
            (x * x + 2.0 * x + 3.0) / (4.0 * x * x + 5.0 * x + 6.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestRational Failed");
        }
    }   
    {
        double x = 0.5, y[1];
        df(x, y);
        auto sq = [](const double x){return x*x;};
        double ref[1] = {
            -(3.0 * x * x + 12.0 * x + 3.0) /
            sq(4.0 * x * x + 5.0 * x + 6.0)
        };
        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestRational Failed");
        }                        
    }
}

void TestTrigonometric() {
    // y = sin(x) + cos(x) + tan(x)
    // dy/dx = cos(x) - sin(x) + sec(x)^2
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto y = std::make_shared<NamedAssignment>("y", 0, sin(x) + cos(x) + tan(x));
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_trigonometric");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {
            std::sin(x) + std::cos(x) + std::tan(x)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestTrigonometric Failed");
        }
    }   
    {
        double x = 0.5, y[1];
        df(x, y);
        auto sq = [](const double x){return x*x;};
        double ref[1] = {
            std::cos(x) - std::sin(x) + (1.0 / sq(std::cos(x)))
        };
        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestTrigonometric Failed");
        }                        
    }       
}

void TestMultiInput() {
    // y = x1 (2x0^2 + x0 + 0.5)
    // dy/dx0 = x1 (4x0 + 1)
    // dy/dx1 = 2x0^2 + x0 + 0.5
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x", 2);
    auto x0 = xArg->GetExpr(0);
    auto x1 = xArg->GetExpr(1);
    auto x0Sq = x0 * x0;
    auto y = std::make_shared<NamedAssignment>("y", 0, x1 * (2.0 * x0Sq + x0 + 0.5));    
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_multiinput");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x0}, {y, x1}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double *, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double *, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x[2] = {0.5, 0.3}, y[1];
        f(x, y);
        double ref[1] = {x[1] * (2.0 * x[0] * x[0] + x[0] + 0.5)};
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestMultiInput Failed");
        }
    }   
    {
        double x[2] = {0.5, 0.3}, y[2];
        df(x, y);
        double ref[2] = {
            x[1] * (4.0 * x[0] + 1.0),
            2.0 * x[0] * x[0] + x[0] + 0.5,
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestMultiInput Failed");
            }                
        }
    }       
}

void TestMultiOutput() {
    // y0 = x^2 + 2x + 3
    // y1 = 4x^2 + 5x + 6
    // dy0/dx = 2x + 2
    // dy1/dx = 8x + 5
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;
    auto y0 = std::make_shared<NamedAssignment>("y", 0, xSq + 2.0 * x + 3.0);
    auto y1 = std::make_shared<NamedAssignment>("y", 1, 4.0 * xSq + 5.0 * x + 6.0);
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y0, y1};
    
    Library lib("func_multioutput");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y0, x}, {y1, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[2];
        f(x, y);
        double ref[2] = {
            x * x + 2.0 * x + 3.0,
            4.0 * x * x + 5.0 * x + 6.0,
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestMultiOutput Failed");
            }
        }
    }   
    {
        double x = 0.5, y[2];
        df(x, y);
        double ref[2] = {
            2.0 * x + 2.0,
            8.0 * x + 5.0
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestMultiOutput Failed");
            }                
        }
    }
}

void TestMultiInputOutput() {
    // y0 = x0^2*x1 + 2x0 + 3x1
    // y1 = 4x0*x1^2 + 5x1^2 + 6x0
    // dy0/dx0 = 2x0*x1 + 2
    // dy1/dx0 = 4x1^2 + 6
    // dy0/dx1 = x0^2 + 3
    // dy1/dx1 = 8*x0*x1 + 10*x1    
    ClearExpressionCache();    
    auto xArg = std::make_shared<DoubleArgument>("x", 2);
    auto x0 = xArg->GetExpr(0);
    auto x1 = xArg->GetExpr(1);
    auto x0Sq = x0 * x0;
    auto x1Sq = x1 * x1;
    auto y0 = std::make_shared<NamedAssignment>("y", 0, x0Sq * x1 + 2.0 * x0 + 3.0 * x1);
    auto y1 = std::make_shared<NamedAssignment>("y", 1, 4.0 * x0 * x1Sq + 5.0 * x1Sq + 6.0 * x0);
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y0, y1};
    
    Library lib("func_multiinout");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y0, x0}, {y1, x0}, {y0, x1}, {y1, x1}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double *, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double *, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x[2] = {0.5, 0.3}, y[2];
        f(x, y);
        double ref[2] = {
            x[0] * x[0] * x[1] + 2.0 * x[0] + 3.0 * x[1],
            4.0 * x[0] * x[1] * x[1] + 5.0 * x[1] * x[1] + 6.0 * x[0],
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestMultiInputOutput Failed");
            }            
        }
    }   
    {
        double x[2] = {0.5, 0.3}, y[4];
        df(x, y);
        double ref[4] = {
            2.0 * x[0] * x[1] + 2.0,
            4.0 * x[1] * x[1] + 6.0,
            x[0] * x[0] + 3.0,
            8.0 * x[0] * x[1] + 10.0 * x[1]
        };
        for (int i = 0; i < 4; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestMultiInputOutput Failed");
            }                
        }
    }       
}

void TestSecondDerivative() {
    // y = 10x^3 + 9x^2 + 8x + 7
    // dy/dx = 30x^2 + 18x + 8
    // d^2y/dx^2 = 60x + 18
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;
    auto xCu = xSq * x;
    auto y = std::make_shared<NamedAssignment>("y", 0, 10.0 * xCu + 9.0 * xSq + 8.0 * x + 7.0);    
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_second_derv");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    std::vector<std::shared_ptr<Expression>> derv2 = Derivatives({{dervOutput[0], x}});
    std::vector<std::shared_ptr<NamedAssignment>> derv2Output;
    for (int i = 0; i < (int)derv2.size(); i++) {
        derv2Output.push_back(std::make_shared<NamedAssignment>("d2ydx2", i, derv2[i]));        
    }
    std::string derv2FuncName = "d2f";    
    lib.AddFunction(input, derv2Output, derv2FuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    
    typedef void (*d2func_t)(const double, double *);
    dfunc_t d2f = (d2func_t)lib.LoadFunction(derv2FuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {
            10.0 * x * x * x + 9.0 * x * x + 8.0 * x + 7.0
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestSecondDerivative Failed");
        }
    }   
    {
        double x = 0.5, y[1];
        df(x, y);
        double ref[1] = {
            30.0 * x * x + 18.0 * x + 8.0
        };
        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestSecondDerivative Failed");
        }                        
    }           
    {
        double x = 0.5, y[1];
        d2f(x, y);
        double ref[1] = {
            60.0 * x + 18.0
        };
        
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestSecondDerivative Failed");
        }                        
    }               
}

void TestJacobian() {
    // f0 = x^2y
    // f0 = 5x + sin(y)
    // J = | 2xy    x^2  |
    //     |  5   cos(y) |
    // det(J) = 2xycos(y) - 5x^2
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto yArg = std::make_shared<DoubleArgument>("y");
    auto x = xArg->GetExpr();    
    auto y = yArg->GetExpr();
    auto xSq = x * x;    
    auto f0 = std::make_shared<NamedAssignment>("f", 0, xSq * y);
    auto f1 = std::make_shared<NamedAssignment>("f", 1, 5.0 * x + sin(y));
    std::vector<std::shared_ptr<Argument>> input = {xArg, yArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {f0, f1};
    
    Library lib("func_jacobian");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{f0, x}, {f0, y}, {f1, x}, {f1, y}});
    auto J = std::make_shared<NamedAssignment>("j", 0, derv[0] * derv[3] - derv[1] * derv[2]);
    std::string jacobianName = "J";    
    lib.AddFunction(input, {J}, jacobianName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, const double, double *);    
    func_t func = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*jac_t)(const double, const double, double *);
    jac_t jac = (jac_t)lib.LoadFunction(jacobianName.c_str());    

    {
        double x = 0.5, y = 0.3, f[2];
        func(x, y, f);
        double ref[2] = {
            x * x * y,
            5.0 * x + std::sin(y)
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(f[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestJacobian Failed");
            }            
        }
    }   
    {
        double x = 0.5, y = 0.3, j[1];
        jac(x, y, j);
        double ref[1] = {
            2.0 * x * y * std::cos(y) - 5.0 * x * x
        };
        for (int i = 0; i < 1; i++) {
            if (fabs(j[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestJacobian Failed");
            }                
        }
    }
}

void TestIfElse() {
    // y = {x^2 + 2x + 3    if x > 0   
    //      4*x^2 + 5x + 6  otherwise}
    // dy/dx = {2x + 2      if x > 0
    //          8x + 5      otherwise}
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto xSq = x * x;    
    auto y = std::make_shared<NamedAssignment>("y", 0, 
        IfElse(Gt(x, 0.0), xSq + 2.0 * x + 3.0, 4.0 * xSq + 5.0 * x + 6.0));
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_ifelse");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, double *);    
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());    
    typedef void (*dfunc_t)(const double, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x = 0.5, y[1];
        f(x, y);
        double ref[1] = {
            x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                      (4.0*x*x + 5.0*x + 6.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElse Failed");
        }
    }
    {
        double x = -0.5, y[1];
        f(x, y);
        double ref[1] = {
            x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                      (4.0*x*x + 5.0*x + 6.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElse Failed");
        }
    }    
    {
        double x = 0.5, y[1];
        df(x, y);
        double ref[1] = {
            x > 0.0 ? (2.0*x + 2.0) : 
                      (8.0*x + 5.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElse Failed");
        }
    }
    {
        double x = -0.5, y[1];
        df(x, y);
        double ref[1] = {
            x > 0.0 ? (2.0*x + 2.0) : 
                      (8.0*x + 5.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElse Failed");
        }
    }
}

void TestIfElseRational() {
    // y = {(x0^2 + 2x0 + 3) / x1  if x1 > 1
    //       x0^2 + 2x0 + 3        otherwise}
    // dy/dx0 = {2x0 + 2 / x1      if x1 > 1
    //           2x0 + 2           otherwise}
    // dy/dx1 = {-(2x0 + 2 / x1^2) if x1 > 1
    //           0                 otherwise}
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x", 2);
    auto x0 = xArg->GetExpr(0);
    auto x1 = xArg->GetExpr(1);
    auto x0Sq = x0 * x0;
    auto y = std::make_shared<NamedAssignment>("y", 0, 
        IfElse(Gt(x1, 1.0), (x0Sq + 2.0 * x0 + 3.0) / x1, (x0Sq + 2.0 * x0 + 3.0)));
    std::vector<std::shared_ptr<Argument>> input = {xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y};
    
    Library lib("func_ifelserational");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);        
    std::vector<std::shared_ptr<Expression>> derv = Derivatives({{y, x0}, {y, x1}});
    std::vector<std::shared_ptr<NamedAssignment>> dervOutput;
    for (int i = 0; i < (int)derv.size(); i++) {
        dervOutput.push_back(std::make_shared<NamedAssignment>("dydx", i, derv[i]));        
    }
    std::string dervFuncName = "df";    
    lib.AddFunction(input, dervOutput, dervFuncName); 
    
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double *, double *);
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());
    typedef void (*dfunc_t)(const double *, double *);
    dfunc_t df = (dfunc_t)lib.LoadFunction(dervFuncName.c_str());    

    {
        double x[2] = {0.5, 2.0}, y[1];
        f(x, y);
        double ref[1] = {
            x[1] > 1.0 ? ((x[0]*x[0] + 2.0*x[0] + 3.0) / x[1]) : 
                          (x[0]*x[0] + 2.0*x[0] + 3.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseRational Failed");
        }
    }
    {
        double x[2] = {0.5, 0.0}, y[1];
        f(x, y);
        double ref[1] = {
            x[1] > 1.0 ? ((x[0]*x[0] + 2.0*x[0] + 3.0) / x[1]) : 
                          (x[0]*x[0] + 2.0*x[0] + 3.0)
        };
        if (fabs(y[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseRational Failed");
        }
    }    
    {
        double x[2] = {0.5, 2.0}, y[2];
        df(x, y);
        double ref[2] = {
            x[1] > 1.0 ? ((2.0*x[0] + 2.0) / x[1]) :
                          (2.0*x[0] + 2.0),
            x[1] > 1.0 ? (-(x[0]*x[0] + 2.0*x[0] + 3.0) / (x[1] * x[1])) :
                         0.0
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestIfElseRational Failed");
            }
        }
    }
    {
        double x[2] = {0.5, 0.0}, y[2];
        df(x, y);
        double ref[2] = {
            x[1] > 1.0 ? ((2.0*x[0] + 2.0) / x[1]) :
                          (2.0*x[0] + 2.0),
            x[1] > 1.0 ? (-(x[0]*x[0] + 2.0*x[0] + 3.0) / (x[1] * x[1])) :
                         0.0
        };
        for (int i = 0; i < 2; i++) {
            if (fabs(y[i] - ref[i]) > 1e-6) {
                throw std::runtime_error("TestIfElseRational Failed");
            }
        }
    }
}

void TestIfElseCycle() {
    // y0 = {x^2 + 2x + 3    if x > 0   
    //       4*x^2 + 5x + 6  otherwise}
    // z  = {5 * y0          if t > 0
    //       x + y0          otherwise}
    // y1 = {z * z           if x > 0
    //       z + z           otherwise}
    // z depends on y0, y1 depends on z
    // y1 must NOT merge with y0
    ClearExpressionCache();
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto tArg = std::make_shared<DoubleArgument>("t");
    auto x = xArg->GetExpr();
    auto t = tArg->GetExpr();
    auto xSq = x * x;
    auto y0 = IfElse(Gt(x, 0.0),       xSq + 2.0 * x + 3.0, 
                                 4.0 * xSq + 5.0 * x + 6.0);
    auto z  = IfElse(Gt(t, 0.0), 5.0 * y0, x + y0);
    auto y1 = std::make_shared<NamedAssignment>("y1", 0, IfElse(Gt(x, 0.0), z * z, z + z));
    std::vector<std::shared_ptr<Argument>> input = {xArg, tArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {y1};
    
    Library lib("func_ifelse_cycle");
    std::string funcName = "f";
    lib.AddFunction(input, output, funcName);
    lib.CompileAndLoad();
    
    typedef void (*func_t)(const double, const double, double *);
    func_t f = (func_t)lib.LoadFunction(funcName.c_str());

    {
        double x = 0.5, t = 0.5, y1[1];
        f(x, t, y1);
        double y0 = x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                             (4.0*x*x + 5.0*x + 6.0);
        double z  = t > 0.0 ? 5.0 * y0 : x + y0;
        double ref[1] = {
            x > 0.0 ? (z * z) : (z + z)
        };
        if (fabs(y1[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseCycle Failed");
        }
    }
    {
        double x = -0.5, t = 0.5, y1[1];
        f(x, t, y1);
        double y0 = x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                             (4.0*x*x + 5.0*x + 6.0);
        double z  = t > 0.0 ? 5.0 * y0 : x + y0;
        double ref[1] = {
            x > 0.0 ? (z * z) : (z + z)
        };
        if (fabs(y1[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseCycle Failed");
        }
    }
    {
        double x = 0.5, t = -0.5, y1[1];
        f(x, t, y1);
        double y0 = x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                             (4.0*x*x + 5.0*x + 6.0);
        double z  = t > 0.0 ? 5.0 * y0 : x + y0;
        double ref[1] = {
            x > 0.0 ? (z * z) : (z + z)
        };
        if (fabs(y1[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseCycle Failed");
        }
    }
    {
        double x = -0.5, t = -0.5, y1[1];
        f(x, t, y1);
        double y0 = x > 0.0 ? (x*x + 2.0*x + 3.0) : 
                             (4.0*x*x + 5.0*x + 6.0);
        double z  = t > 0.0 ? 5.0 * y0 : x + y0;
        double ref[1] = {
            x > 0.0 ? (z * z) : (z + z)
        };
        if (fabs(y1[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestIfElseCycle Failed");
        }
    }
}

void TestStruct() {
    // struct Foo {
    //     double x;
    //     double y[2];
    // } foo;
    // z = foo.x + foo.y[0] * foo.y[1];
    ClearExpressionCache();
    StructType fooType("Foo", {
        std::make_shared<DoubleArgument>("x"),
        std::make_shared<DoubleArgument>("y", 2)});
    auto fooArg = std::make_shared<StructArgument>("foo", fooType);
    auto x = fooArg->GetArg("x")->GetExpr();
    auto y0 = fooArg->GetArg("y")->GetExpr(0);
    auto y1 = fooArg->GetArg("y")->GetExpr(1);
    auto z = std::make_shared<NamedAssignment>("z", 0, x + y0 * y1);
    std::vector<std::shared_ptr<Argument>> input = {fooArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {z};
    
    Library lib("func_struct");
    lib.AddStruct(fooType);
    lib.AddFunction(input, output, "f");
    lib.CompileAndLoad();
    // Is it possible/worthwhile to use some template/macro magic to generate StructType from this?
    struct Foo {
        double x;
        double y[2];
    };
    typedef void (*func_t)(const Foo *, double *);
    func_t f = (func_t)lib.LoadFunction("f");
    {
        Foo foo = {0.5, {0.6, 0.7}};
        double z[1];
        f(&foo, z);
        double ref[1] = {
            foo.x + foo.y[0] * foo.y[1]
        };
        if (fabs(z[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestStruct Failed");
        }
    }
}

void TestUnion() {
    // struct Foo {
    //     double a;
    //     double b[2];
    // } foo;
    // struct Bar {
    //     double c[2];
    //     double d;
    // } bar;
    // union FooBar {
    //     Foo foo;
    //     Bar bar;
    // } foobar;
    // z = x > 0 ? foo.a * foo.b[0] * foo.b[1] : bar.c[0] + bar.c[1] + bar.d;
    ClearExpressionCache();
    StructType fooType("Foo", {
        std::make_shared<DoubleArgument>("a"),
        std::make_shared<DoubleArgument>("b", 2)});
    StructType barType("Bar", {
        std::make_shared<DoubleArgument>("c", 2),
        std::make_shared<DoubleArgument>("d")});
    StructType foobarType("FooBar", {
        std::make_shared<StructArgument>("foo", fooType),
        std::make_shared<StructArgument>("bar", barType)},
        true);
    auto fooBarArg = std::make_shared<StructArgument>("foobar", foobarType);
    auto xArg = std::make_shared<DoubleArgument>("x");
    auto x = xArg->GetExpr();
    auto foo = fooBarArg->GetArg("foo")->GetArg("a")->GetExpr() *
               fooBarArg->GetArg("foo")->GetArg("b")->GetExpr(0) *
               fooBarArg->GetArg("foo")->GetArg("b")->GetExpr(1);
    auto bar = fooBarArg->GetArg("bar")->GetArg("c")->GetExpr(0) +
               fooBarArg->GetArg("bar")->GetArg("c")->GetExpr(1) +
               fooBarArg->GetArg("bar")->GetArg("d")->GetExpr();
    auto z = std::make_shared<NamedAssignment>("z", 0, IfElse(Gt(x, 0.0), foo, bar));
    std::vector<std::shared_ptr<Argument>> input = {fooBarArg, xArg};
    std::vector<std::shared_ptr<NamedAssignment>> output = {z};
    
    Library lib("func_union");
    lib.AddStruct(fooType);
    lib.AddStruct(barType);
    lib.AddStruct(foobarType);
    lib.AddFunction(input, output, "f");
    lib.CompileAndLoad();
    struct Foo {
        double a;
        double b[2];
    };
    struct Bar {
        double c[2];
        double d;
    };
    union FooBar {
        Foo foo;
        Bar bar;
    };
    typedef void (*func_t)(const FooBar *, const double, double *);
    func_t f = (func_t)lib.LoadFunction("f");
    {
        double x = -0.3;
        FooBar foobar;
        foobar.bar = {{0.6, 0.7}, 0.8};
        double z[1];
        f(&foobar, x, z);
        double ref[1] = {
            x > 0.0 ? foobar.foo.a * foobar.foo.b[0] * foobar.foo.b[1] :
                      foobar.bar.c[0] + foobar.bar.c[1] + foobar.bar.d
        };
        if (fabs(z[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestUnion Failed");
        }
    }
    {
        double x = 0.3;
        FooBar foobar;
        foobar.foo = {0.3, {0.4, 0.5}};
        double z[1];
        f(&foobar, x, z);
        double ref[1] = {
            x > 0.0 ? foobar.foo.a * foobar.foo.b[0] * foobar.foo.b[1] :
                      foobar.bar.c[0] + foobar.bar.c[1] + foobar.bar.d
        };
        if (fabs(z[0] - ref[0]) > 1e-6) {
            throw std::runtime_error("TestUnion Failed");
        }
    }
}

int main(int argc, char *argv[]) {
    try {
        TestConstant();
        std::cerr << "TestConstant Passed" << std::endl;
        TestLinear();
        std::cerr << "TestLinear Passed" << std::endl;
        TestPolynomial();
        std::cerr << "TestPolynomial Passed" << std::endl;        
        TestRational();
        std::cerr << "TestRational Passed" << std::endl;
        TestTrigonometric();
        std::cerr << "TestTrigonometric Passed" << std::endl;
        TestMultiInput();        
        std::cerr << "TestMultiInput Passed" << std::endl;
        TestMultiOutput();
        std::cerr << "TestMultiOutput Passed" << std::endl;
        TestMultiInputOutput();
        std::cerr << "TestMultiInputOutput Passed" << std::endl;
        TestSecondDerivative();
        std::cerr << "TestSecondDerivative Passed" << std::endl;    
        TestJacobian();
        std::cerr << "TestJacobian Passed" << std::endl;
        TestIfElse();
        std::cerr << "TestIfElse Passed" << std::endl;
        TestIfElseRational();
        std::cerr << "TestIfElseRational Passed" << std::endl;
        TestIfElseCycle();
        std::cerr << "TestIfElseCycle Passed" << std::endl;
        TestStruct();
        std::cerr << "TestStruct Passed" << std::endl;
        TestUnion();
        std::cerr << "TestUnion Passed" << std::endl;
    } catch (std::exception &ex) {
        std::cout << ex.what() << std::endl;
        return -1;
    }
    std::cout << "All Passed" << std::endl;
    return 0;
}
