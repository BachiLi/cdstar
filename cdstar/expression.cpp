#include "expression.h"
#include "derivativegraph.h"
#include "argument.h"

#include <iostream>

namespace cdstar {

struct ExprHasher {
    std::size_t operator()(const std::shared_ptr<Expression> &expr) const {
        return expr->GetHash();
    }
};

struct ExprComparator {
    bool operator()(const std::shared_ptr<Expression> &expr0,
                    const std::shared_ptr<Expression> &expr1) const {
        if (expr0 == expr1) {
            return true;
        }
        return expr0->Equal(expr1);
    }
};

std::unordered_set<std::shared_ptr<Expression>, ExprHasher, ExprComparator> g_ExprCache;

std::ostream& PrintTab(const int tabNum, std::ostream &os) {
    for (int i = 0; i < tabNum; i++)
        os << "\t";
    return os;
}

////////////////////////////////////////////////////////////////////
void Variable::Print() const {
    if (m_Index >= 0) {
        std::cerr << m_Name << "[" << m_Index << "]";
    } else {
        std::cerr << m_Name << std::endl;
    }
}

void Variable::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
}

std::string Variable::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    if (m_Index >= 0) {
        return m_Name + std::string("[") + std::to_string(m_Index) + std::string("]");
    } else {
        return m_Name;
    }
}

std::vector<std::shared_ptr<Expression>> Variable::Children() const {
    return {};
}

std::vector<std::shared_ptr<Expression>> Variable::Dervs() const {
    return {};
}

void Constant::Print() const {    
    std::cerr << m_Value;
}

void Constant::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
}

std::string Constant::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::to_string(m_Value);
}

std::vector<std::shared_ptr<Expression>> Constant::Children() const {
    return {};
}

std::vector<std::shared_ptr<Expression>> Constant::Dervs() const {
    return {};
}
////////////////////////////////////////////////////////////////////////////
void IntegerConstant::Print() const {    
    std::cerr << m_Value;
}

void IntegerConstant::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
}

std::string IntegerConstant::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::to_string(m_Value);
}

std::vector<std::shared_ptr<Expression>> IntegerConstant::Children() const {
    return {};
}

std::vector<std::shared_ptr<Expression>> IntegerConstant::Dervs() const {
    return {};
}
////////////////////////////////////////////////////////////////////////////
void NamedAssignment::Print() const {
    std::cerr << m_Name << "[" << m_Index << "] = ";
    m_Expr->Print();
}

void NamedAssignment::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = " << m_Expr->GetEmitName(exprMap) << ";" << std::endl;
}

std::string NamedAssignment::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return m_Name + std::string("[") + std::to_string(m_Index) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> NamedAssignment::Children() const {
    return {m_Expr};
}

std::vector<std::shared_ptr<Expression>> NamedAssignment::Dervs() const {
    return {std::make_shared<Constant>(1.0)};
}
////////////////////////////////////////////////////////////////////////////

std::string UnaryAssignment::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::string("_t[") + std::to_string(exprMap.at(this).index) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> UnaryAssignment::Children() const {
    return {m_Expr};
}
////////////////////////////////////////////////////////////////////////////
void Negate::Print() const {
    std::cerr << "(-";
    m_Expr->Print();
    std::cerr << ")";
}

void Negate::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = -" << m_Expr->GetEmitName(exprMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Negate::Dervs() const {
    return {std::make_shared<Constant>(-1.0)};
}
////////////////////////////////////////////////////////////////////////////
void Inverse::Print() const {
    std::cerr << "(1.0 / ";
    m_Expr->Print();
    std::cerr << ")";
}

void Inverse::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = 1.0 / " << m_Expr->GetEmitName(exprMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Inverse::Dervs() const {
    // d(1/x)/dx = -1/x^2
    std::shared_ptr<Inverse> thisPtr = std::const_pointer_cast<Inverse>(shared_from_this());
    return {-(thisPtr * thisPtr)};
}
////////////////////////////////////////////////////////////////////////////
void Sin::Print() const {
    std::cerr << "sin(";
    m_Expr->Print();
    std::cerr << ")";
}

void Sin::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = sin(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Sin::Dervs() const {
    return {cos(m_Expr)};
}
////////////////////////////////////////////////////////////////////////////
void Cos::Print() const {
    std::cerr << "cos(";
    m_Expr->Print();
    std::cerr << ")";
}

void Cos::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = cos(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Cos::Dervs() const {
    return {-sin(m_Expr)};
}
////////////////////////////////////////////////////////////////////////////
void Tan::Print() const {
    std::cerr << "tan(";
    m_Expr->Print();
    std::cerr << ")";
}

void Tan::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = tan(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Tan::Dervs() const {
    auto c = cos(m_Expr);
    return {Inv(c * c)};
}
////////////////////////////////////////////////////////////////////////////
void Sqrt::Print() const {
    std::cerr << "sqrt(";
    m_Expr->Print();
    std::cerr << ")";
}

void Sqrt::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = sqrt(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Sqrt::Dervs() const {
    std::shared_ptr<Sqrt> thisPtr = std::const_pointer_cast<Sqrt>(shared_from_this());
    return {Inv(2.0 * thisPtr)};
}
////////////////////////////////////////////////////////////////////////////
void ASin::Print() const {
    std::cerr << "asin(";
    m_Expr->Print();
    std::cerr << ")";
}

void ASin::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = asin(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> ASin::Dervs() const {
    return {Inv(sqrt(1.0 - m_Expr * m_Expr))};
}
////////////////////////////////////////////////////////////////////////////
void ACos::Print() const {
    std::cerr << "acos(";
    m_Expr->Print();
    std::cerr << ")";
}

void ACos::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = acos(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> ACos::Dervs() const {
    return {-Inv(sqrt(1.0 - m_Expr * m_Expr))};
}
////////////////////////////////////////////////////////////////////////////
void Log::Print() const {
    std::cerr << "log(";
    m_Expr->Print();
    std::cerr << ")";
}

void Log::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = log(" << m_Expr->GetEmitName(exprMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Log::Dervs() const {
    return {Inv(m_Expr)};
}
////////////////////////////////////////////////////////////////////////////
std::string BinaryAssignment::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::string("_t[") + std::to_string(exprMap.at(this).index) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> BinaryAssignment::Children() const {
    return {m_Expr0, m_Expr1};
}
////////////////////////////////////////////////////////////////////////////
void Add::Print() const {
    std::cerr << "(";
    m_Expr0->Print();
    std::cerr << " + ";
    m_Expr1->Print();
    std::cerr << ")";
}

void Add::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = " << 
        m_Expr0->GetEmitName(exprMap) << " + " <<
        m_Expr1->GetEmitName(exprMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Add::Dervs() const {
    auto c = std::make_shared<Constant>(1.0);
    return {c, c};
}
////////////////////////////////////////////////////////////////////////////
void Multiply::Print() const {
    m_Expr0->Print();
    std::cerr << " * ";
    m_Expr1->Print();
}

void Multiply::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << GetEmitName(exprMap) << " = " << 
        m_Expr0->GetEmitName(exprMap) << " * " <<
        m_Expr1->GetEmitName(exprMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Multiply::Dervs() const {
    return {m_Expr1, m_Expr0};
}    
////////////////////////////////////////////////////////////////////////////
std::string Boolean::OpToString() const {
    switch(m_Op) {
        case GREATER:
            return " > ";
            break;
        case GREATER_OR_EQUAL:
            return " >= ";
            break;
        case EQUAL:
            return " == ";
            break;
        case NOT_EQUAL:
            return " != ";
            break;
        case LESS_OR_EQUAL:
            return " <= ";
            break;
        case LESS:
            return " < ";
            break;
    }    
    return "";
}

Boolean::Op Boolean::ReverseOp() const {
    switch(m_Op) {
        case GREATER:
            return LESS_OR_EQUAL;
            break;
        case GREATER_OR_EQUAL:
            return LESS;
            break;
        case EQUAL:
            return NOT_EQUAL;
            break;
        case NOT_EQUAL:
            return EQUAL;
            break;
        case LESS_OR_EQUAL:
            return LESS;
            break;
        case LESS:
            return LESS_OR_EQUAL;
            break;
    }    
    return EQUAL;    
}

void Boolean::Print() const {
    std::cerr << "(";
    m_Expr0->Print();
    std::cerr << OpToString();
    m_Expr1->Print();
    std::cerr << ")";
}

void Boolean::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    os << "int " << GetEmitName(exprMap) << " = " <<
        m_Expr0->GetEmitName(exprMap) << OpToString() << m_Expr1->GetEmitName(exprMap) << ";" << std::endl;
}

std::string Boolean::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::string("b") + std::to_string(exprMap.at(this).index);
}

std::vector<std::shared_ptr<Expression>> Boolean::Children() const {
    return {m_Expr0, m_Expr1};
}

std::vector<std::shared_ptr<Expression>> Boolean::Dervs() const {
    auto zero = std::make_shared<Constant>(0.0);
    return {zero, zero};
}
////////////////////////////////////////////////////////////////////////////
void CondExpr::Print() const {
    std::cerr << "(";
    m_Cond->Print();
    std::cerr << " ? ";
    m_TrueExpr->Print();
    std::cerr << " : ";
    m_FalseExpr->Print();
    std::cerr << ")";
}

void CondExpr::Emit(const CFGBlock *block, const std::unordered_map<const Expression*, ExprInfo> &exprMap, std::ostream &os) const {
    if (block->condId == 0) { // emit true
        os << GetEmitName(exprMap) << " = " << m_TrueExpr->GetEmitName(exprMap) << ";" << std::endl;
    } else { // emit false
        os << GetEmitName(exprMap) << " = " << m_FalseExpr->GetEmitName(exprMap) << ";" << std::endl;
    }
}

std::string CondExpr::GetEmitName(const std::unordered_map<const Expression*, ExprInfo> &exprMap) const {
    return std::string("_t[") + std::to_string(exprMap.at(this).index) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> CondExpr::Children() const {
    return {m_Cond, m_TrueExpr, m_FalseExpr};
}

std::vector<std::shared_ptr<Expression>> CondExpr::Dervs() const {
    auto one  = std::make_shared<Constant>(1.0);
    auto zero = std::make_shared<Constant>(0.0);
    return {zero,
            IfElse(m_Cond, one, zero),
            IfElse(m_Cond, zero, one)};
}

std::shared_ptr<Expression> operator+(const std::shared_ptr<Expression> expr0, 
                                      const std::shared_ptr<Expression> expr1) {
    if (expr0->Type() == ET_CONSTANT) {
        return expr0->GetConstant() + expr1;
    }
    if (expr1->Type() == ET_CONSTANT) {
        return expr0 + expr1->GetConstant();
    }
    if (expr0->Type() == ET_MULTIPLY && expr1->Type() == ET_MULTIPLY) {
        auto children0 = expr0->Children();
        auto children1 = expr1->Children();
        if (children0[0] == children1[0]) {
            return children0[0] * (children0[1] + children1[1]);
        } else if (children0[0] == children1[1]) {
            return children0[0] * (children0[1] + children1[0]);
        } else if (children0[1] == children1[0]) {
            return children0[1] * (children0[0] + children1[1]);
        } else if (children0[1] == children1[1]) {
            return children0[1] * (children0[0] + children1[0]);
        }
    }

    if (expr0->Type() == ET_CONDEXPR && expr1->Type() == ET_CONDEXPR) {
        auto children0 = expr0->Children();
        auto children1 = expr1->Children();
        if (children0[0] == children1[0]) {
            return IfElse(std::dynamic_pointer_cast<Boolean>(children0[0]),
                children0[1] + children1[1], children0[2] + children1[2]);
        }
    }

    if (expr0->Type() == ET_CONDEXPR) {
        auto children = expr0->Children();
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), expr1, children[2] + expr1);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), children[1] + expr1, expr1);
        }
    }
    if (expr1->Type() == ET_CONDEXPR) {
        auto children = expr1->Children();
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), expr0, expr0 + children[2]);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), expr0 + children[1], expr0);
        }
    }
    auto ret0 = std::make_shared<Add>(expr0, expr1);
    auto ret1 = std::make_shared<Add>(expr1, expr0);
    return CacheExpression(ret0, ret1);
}

std::shared_ptr<Expression> operator+(const double v0,
                                      const std::shared_ptr<Expression> expr1) {
    if (v0 == 0.0) {
        return expr1;
    }
    if (expr1->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(v0 + expr1->GetConstant());
    }
    if (expr1->Type() == ET_ADD && expr1->HasPartialConstant()) {
        auto children = expr1->Children();
        if (expr1->HasPartialConstant() == 1) {
            return (v0 + children[0]->GetConstant()) + children[1];
        } else {
            return (v0 + children[1]->GetConstant()) + children[0];
        }
    }
    if (expr1->Type() == ET_CONDEXPR) {
        auto children = expr1->Children();
        if (children[1]->Type() == ET_CONSTANT) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(),
                          std::make_shared<Constant>(v0 + children[1]->GetConstant()),
                          v0 + children[2]);
        }
        if (children[2]->Type() == ET_CONSTANT) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(),
                          v0 + children[1], 
                          std::make_shared<Constant>(v0 + children[2]->GetConstant()));
        }
    }
    auto c = std::make_shared<Constant>(v0);
    auto ret0 = std::make_shared<Add>(c, expr1);
    auto ret1 = std::make_shared<Add>(expr1, c);
    return CacheExpression(ret0, ret1);
}

std::shared_ptr<Expression> operator-(const std::shared_ptr<Expression> expr) {
    if (expr->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(-expr->GetConstant());
    }
    if (expr->Type() == ET_NEGATE) {
        return expr->Children()[0];
    }
    return CacheExpression(std::make_shared<Negate>(expr));
}

std::shared_ptr<Expression> operator-(const std::shared_ptr<Expression> expr0,
                                      const std::shared_ptr<Expression> expr1) {
    if (expr0.get() == expr1.get()) {
        return std::make_shared<Constant>(0.0);
    }
    return expr0 + (-expr1);
}

std::shared_ptr<Expression> operator*(const std::shared_ptr<Expression> expr0, 
                                      const std::shared_ptr<Expression> expr1) {
    if (expr0->Type() == ET_CONSTANT) {
        return expr0->GetConstant() * expr1;
    }
    if (expr1->Type() == ET_CONSTANT) {
        return expr0 * expr1->GetConstant();
    }

    if (expr0->Type() == ET_CONDEXPR) {
        auto children = expr0->Children();
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), 0.0, children[2] * expr1);
        }
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 1.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), expr1, children[2] * expr1);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), children[1] * expr1, 0.0);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 1.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr0);
            return IfElse(condExpr->GetCond(), children[1] * expr1, expr1);
        }
    }
    if (expr1->Type() == ET_CONDEXPR) {
        auto children = expr1->Children();
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), 0.0, expr0 * children[2]);
        }
        if (children[1]->Type() == ET_CONSTANT && children[1]->GetConstant() == 1.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), expr0, expr0 * children[2]);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 0.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), expr0 * children[1], 0.0);
        }
        if (children[2]->Type() == ET_CONSTANT && children[2]->GetConstant() == 1.0) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(), expr0 * children[1], expr0);
        }
    }
    auto ret0 = std::make_shared<Multiply>(expr0, expr1);
    auto ret1 = std::make_shared<Multiply>(expr1, expr0);
    return CacheExpression(ret0, ret1);
}

std::shared_ptr<Expression> operator*(const double v0,
                                      const std::shared_ptr<Expression> expr1) {
    if (v0 == 0.0) {
        return std::make_shared<Constant>(0.0);
    } else if (v0 == 1.0) {
        return expr1;
    } else if (v0 == -1.0) {
        return -expr1;
    }
    if (expr1->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(v0 * expr1->GetConstant());
    }
    if (expr1->Type() == ET_MULTIPLY && expr1->HasPartialConstant()) {
        auto children = expr1->Children();
        if (expr1->HasPartialConstant() == 1) {
            return (v0 * children[0]->GetConstant()) * children[1];
        } else {
            return (v0 * children[1]->GetConstant()) * children[0];
        }
    }

    if (expr1->Type() == ET_CONDEXPR) {
        auto children = expr1->Children();
        if (children[1]->Type() == ET_CONSTANT) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(),
                          v0 * children[1]->GetConstant(),
                          v0 * children[2]);
        }
        if (children[2]->Type() == ET_CONSTANT) {
            auto condExpr = std::dynamic_pointer_cast<CondExpr>(expr1);
            return IfElse(condExpr->GetCond(),
                          v0 * children[1],
                          v0 * children[2]->GetConstant());
        }
    }

    auto c = std::make_shared<Constant>(v0);
    auto ret0 = std::make_shared<Multiply>(c, expr1);
    auto ret1 = std::make_shared<Multiply>(expr1, c);
    return CacheExpression(ret0, ret1);
}

inline std::shared_ptr<Expression> Inv(const std::shared_ptr<Expression> expr) {
    if (expr->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(1.0 / expr->GetConstant());
    }
    if (expr->Type() == ET_INVERSE) {
        return expr->Children()[0];
    }
    return CacheExpression(std::make_shared<Inverse>(expr));
}

std::shared_ptr<Expression> operator/(const std::shared_ptr<Expression> expr0,
                                      const std::shared_ptr<Expression> expr1) {
    if (expr0.get() == expr1.get()) {
        return std::make_shared<Constant>(1.0);
    }
    return expr0 * Inv(expr1);
}

std::shared_ptr<Boolean> Gt(const std::shared_ptr<Expression> expr0,
                            const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::GREATER, expr0, expr1);
    auto ret1 = ret0->Swap();
    // Should fix this?
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Boolean> Gte(const std::shared_ptr<Expression> expr0,
                             const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::GREATER_OR_EQUAL, expr0, expr1);
    auto ret1 = ret0->Swap();
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Boolean> Eq(const std::shared_ptr<Expression> expr0,
                            const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::EQUAL, expr0, expr1);
    auto ret1 = ret0->Swap();
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Boolean> Neq(const std::shared_ptr<Expression> expr0,
                             const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::NOT_EQUAL, expr0, expr1);
    auto ret1 = ret0->Swap();
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Boolean> Lte(const std::shared_ptr<Expression> expr0,
                             const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::LESS_OR_EQUAL, expr0, expr1);
    auto ret1 = ret0->Swap();
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Boolean> Lt(const std::shared_ptr<Expression> expr0,
                            const std::shared_ptr<Expression> expr1) {
    auto ret0 = std::make_shared<Boolean>(Boolean::LESS, expr0, expr1);
    auto ret1 = ret0->Swap();
    return std::dynamic_pointer_cast<Boolean>(CacheExpression(ret0, ret1));
}

std::shared_ptr<Expression> IfElse(const std::shared_ptr<Boolean> cond,
                                   const std::shared_ptr<Expression> trueExpr,
                                   const std::shared_ptr<Expression> falseExpr) {
    if (trueExpr->Type() == ET_CONDEXPR) {
        auto condExpr = std::dynamic_pointer_cast<CondExpr>(trueExpr);
        if (condExpr->GetCond() == cond) {
            return IfElse(cond, condExpr->Children()[1], falseExpr);
        }
    }
    if (falseExpr->Type() == ET_CONDEXPR) {
        auto condExpr = std::dynamic_pointer_cast<CondExpr>(falseExpr);
        if (condExpr->GetCond() == cond) {
            return IfElse(cond, trueExpr, condExpr->Children()[2]);
        }
    }
    auto ret0 = std::make_shared<CondExpr>(cond, trueExpr, falseExpr);
    auto ret1 = ret0->Swap();
    return CacheExpression(ret0, ret1);
}

void PrintExpressionCache() {
    for (auto it : g_ExprCache) {
        it->Print();
        std::cerr << std::endl;
        std::cerr << "use_count:" << it.use_count() << std::endl;
    }
}

void ClearExpressionCache() {
    std::unordered_set<std::shared_ptr<Expression>, ExprHasher, ExprComparator>().swap(g_ExprCache);
}

std::shared_ptr<Expression> CacheExpression(const std::shared_ptr<Expression> &expr) {
    auto it = g_ExprCache.find(expr);
    if (it == g_ExprCache.end()) {
        g_ExprCache.insert(expr);
        return expr;
    }
    return *it;
}

std::shared_ptr<Expression> CacheExpression(const std::shared_ptr<Expression> &expr0,
                                            const std::shared_ptr<Expression> &expr1) {
    auto it = g_ExprCache.find(expr0);
    if (it == g_ExprCache.end()) {
        it = g_ExprCache.find(expr1);
        if (it == g_ExprCache.end()) {
            g_ExprCache.insert(expr0);
            return expr0;
        }
    }

    return *it;
}

std::vector<std::shared_ptr<Expression>> Derivatives(const std::vector<ExprPtrPair> &dervExprs) {
    DerivativeGraph dervGraph(dervExprs);
    return dervGraph.Derivatives();
}

void EmitFunction(const std::vector<std::shared_ptr<Argument>> &inputs,
                  const std::vector<std::shared_ptr<NamedAssignment>> &outputs,
                  const std::string &name,
                  std::ostream &os) {
    std::unordered_map<std::string, int> varMap;
    for (auto &var : outputs) {
        varMap[var->GetName()] = std::max(varMap[var->GetName()], var->GetIndex() + 1);
    }

    os << "void " << name << "(";
    std::unordered_set<std::string> varSet;
    bool first = true;
    for (auto &arg : inputs) {
        if (!first) {
            os << ", ";
        }
        first = false;
        os << "const " << arg->GetDeclaration();
    }
    for (auto &var : outputs) {
        if (varSet.find(var->GetName()) == varSet.end()) {
            if (!first) {
                os << ", ";
            }
            first = false;
            os << "double " << var->GetName() << "[" << varMap[var->GetName()] << "]";
            varSet.insert(var->GetName());
        }
    }
    os << ") {" << std::endl;
    std::shared_ptr<CFGBlock> block = BuildCFGBlock(outputs);
    Optimize(block);
    Emit(block, os);
    os << "}" << std::endl;
}

bool DetectCycle(const std::shared_ptr<Expression> expr, 
                 std::unordered_set<const Expression*> &dfsSet,
                 std::unordered_set<const Expression*> &recSet) {
    if (dfsSet.find(expr.get()) == dfsSet.end()) {
        dfsSet.insert(expr.get());
        recSet.insert(expr.get());

        for (auto child : expr->Children()) {
            if (dfsSet.find(child.get()) != dfsSet.end()) {
                continue;
            }
            bool ret = DetectCycle(child, dfsSet, recSet);
            if (ret) {
                return true;
            }
            if (recSet.find(child.get()) != recSet.end()) {
                return true;
            }
        }
    }
    recSet.erase(expr.get());
    return false;
}

bool DetectCycle(const std::shared_ptr<Expression> expr) {
    std::unordered_set<const Expression*> dfsSet, recSet;
    return DetectCycle(expr, dfsSet, recSet);
}

} //namespace cdstar