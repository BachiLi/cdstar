#include "expression.h"
#include "derivativegraph.h"

#include <iostream>
#include <unordered_set>

namespace cdstar {

struct ExprHasher {
    std::size_t operator()(const std::shared_ptr<Expression> &expr) const {
        return expr->GetHash();
    }
};

struct ExprComparator {
    bool operator()(const std::shared_ptr<Expression> &expr0, 
                    const std::shared_ptr<Expression> &expr1) const {
        return expr0->Equal(expr1);
    }
};

std::unordered_set<std::shared_ptr<Expression>, ExprHasher, ExprComparator> g_ExprCache;

void AssignmentMap::Register(const Expression *expr) {
    if (m_ExprMap.find(expr) == m_ExprMap.end()) {
        m_ExprMap[expr] = m_AssignCount++;            
    }
}

int AssignmentMap::GetIndex(const Expression *expr) const {
    auto it = m_ExprMap.find(expr);
    if (it != m_ExprMap.end()) {
        return it->second;
    }
    return -1;
}

void Variable::Print() const {
    std::cerr << m_Name << "[" << m_Index << "]";
}

void Variable::Register(AssignmentMap &assignMap) const {
}

void Variable::Emit(AssignmentMap &assignMap, std::ostream &os) const {
}

std::string Variable::GetEmitName(const AssignmentMap &assignMap) const {
    return m_Name + std::string("[") + std::to_string(m_Index) + std::string("]");    
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

void Constant::Register(AssignmentMap &assignMap) const {    
}

void Constant::Emit(AssignmentMap &assignMap, std::ostream &os) const {
}

std::string Constant::GetEmitName(const AssignmentMap &assignMap) const {
    return std::to_string(m_Value);    
}

std::vector<std::shared_ptr<Expression>> Constant::Children() const {
    return std::vector<std::shared_ptr<Expression>>();
}

std::vector<std::shared_ptr<Expression>> Constant::Dervs() const {
    return std::vector<std::shared_ptr<Expression>>();
}

void NamedAssignment::Print() const {
    std::cerr << m_Name << "[" << m_Index << "] = ";
    m_Expr->Print();
}

void NamedAssignment::Register(AssignmentMap &assignMap) const {
    m_Expr->Register(assignMap);    
}

void NamedAssignment::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = " << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
    assignMap.SetEmitted(this);
}

std::string NamedAssignment::GetEmitName(const AssignmentMap &assignMap) const {
    return m_Name + std::string("[") + std::to_string(m_Index) + std::string("]");    
}

std::vector<std::shared_ptr<Expression>> NamedAssignment::Children() const {
    return {m_Expr};
}

std::vector<std::shared_ptr<Expression>> NamedAssignment::Dervs() const {
    return {std::make_shared<Constant>(1.0)};
}

void UnaryAssignment::Register(AssignmentMap &assignMap) const {
    m_Expr->Register(assignMap);
    assignMap.Register(this);
}

std::string UnaryAssignment::GetEmitName(const AssignmentMap &assignMap) const {
    return std::string("t[") + std::to_string(assignMap.GetIndex(this)) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> UnaryAssignment::Children() const {
    return {m_Expr};
}

void Negate::Print() const {
    std::cerr << "(-";
    m_Expr->Print();
    std::cerr << ")";
}

void Negate::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = -" << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Negate::Dervs() const {
    return {std::make_shared<Constant>(-1.0)};
}

void Inverse::Print() const {
    std::cerr << "(1.0 / ";
    m_Expr->Print();
    std::cerr << ")";
}

void Inverse::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = 1.0 / " << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Inverse::Dervs() const {
    // d(1/x)/dx = -1/x^2
    return {-(std::shared_ptr<Expression>(const_cast<Inverse*>(this)) * 
              std::shared_ptr<Expression>(const_cast<Inverse*>(this)))};
}

void Sin::Print() const {
    std::cerr << "sin(";
    m_Expr->Print();
    std::cerr << ")";
}

void Sin::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = sin(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Sin::Dervs() const {    
    return {cos(m_Expr)};
}

void Cos::Print() const {
    std::cerr << "cos(";
    m_Expr->Print();
    std::cerr << ")";
}

void Cos::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = cos(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Cos::Dervs() const {    
    return {-sin(m_Expr)};
}

void Tan::Print() const {
    std::cerr << "tan(";
    m_Expr->Print();
    std::cerr << ")";
}

void Tan::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Expr->Emit(assignMap, os);
    os << "\t" << GetEmitName(assignMap) << " = tan(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Tan::Dervs() const {    
    auto c = cos(m_Expr);
    return {Inv(c * c)};
}

void BinaryAssignment::Register(AssignmentMap &assignMap) const {
    m_Expr0->Register(assignMap);
    m_Expr1->Register(assignMap);
    assignMap.Register(this);
}

std::string BinaryAssignment::GetEmitName(const AssignmentMap &assignMap) const {
    return std::string("t[") + std::to_string(assignMap.GetIndex(this)) + std::string("]");
}

std::vector<std::shared_ptr<Expression>> BinaryAssignment::Children() const {    
    std::vector<std::shared_ptr<Expression>> ret;
    ret.push_back(m_Expr0);
    ret.push_back(m_Expr1);
    return ret;
}

void Add::Print() const {
    std::cerr << "(";
    m_Expr0->Print();
    std::cerr << " + ";
    m_Expr1->Print();
    std::cerr << ")";
}

void Add::Emit(AssignmentMap &assignMap, std::ostream &os) const {    
    if (assignMap.IsEmitted(this)) {
        return;
    }    
    m_Expr0->Emit(assignMap, os);
    m_Expr1->Emit(assignMap, os);       
    os << "\t" << GetEmitName(assignMap) << " = " << m_Expr0->GetEmitName(assignMap) << " + " << 
                                                     m_Expr1->GetEmitName(assignMap) << ";" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Add::Dervs() const {    
    auto c = std::make_shared<Constant>(1.0);
    return {c, c};
}

void Multiply::Print() const {
    m_Expr0->Print();
    std::cerr << " * ";
    m_Expr1->Print();
}

void Multiply::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }    
    m_Expr0->Emit(assignMap, os);
    m_Expr1->Emit(assignMap, os);       
    os << "\t" << GetEmitName(assignMap) << " = " << m_Expr0->GetEmitName(assignMap) << " * " << 
                                                     m_Expr1->GetEmitName(assignMap) << ";" << std::endl;
    assignMap.SetEmitted(this);
}

std::vector<std::shared_ptr<Expression>> Multiply::Dervs() const {    
    return {m_Expr1, m_Expr0};
}    

std::shared_ptr<Expression> operator+(const std::shared_ptr<Expression> expr0, 
                                      const std::shared_ptr<Expression> expr1) {
    if (expr0->Type() == ET_CONSTANT && expr1->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(expr0->GetConstant() + expr1->GetConstant());
    }
    if (expr0->Type() == ET_CONSTANT && expr0->GetConstant() == 0.0) {
        return expr1;
    }
    if (expr1->Type() == ET_CONSTANT && expr1->GetConstant() == 0.0) {
        return expr0;
    }
    if (expr0->Type() == ET_ADD && expr0->HasPartialConstant() && expr1->Type() == ET_CONSTANT) {
        auto children = expr0->Children();
        if (expr0->HasPartialConstant() == 1) {
            return (expr1->GetConstant() + children[0]->GetConstant()) + children[1];
        } else {
            return (expr1->GetConstant() + children[1]->GetConstant()) + children[0];
        }
    }
    if (expr1->Type() == ET_ADD && expr1->HasPartialConstant() && expr0->Type() == ET_CONSTANT) {
        auto children = expr1->Children();
        if (expr1->HasPartialConstant() == 1) {
            return (expr0->GetConstant() + children[0]->GetConstant()) + children[1];
        } else {
            return (expr0->GetConstant() + children[1]->GetConstant()) + children[0];
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
    if (expr0->Type() == ET_CONSTANT && expr1->Type() == ET_CONSTANT) {
        return std::make_shared<Constant>(expr0->GetConstant() * expr1->GetConstant());
    }
    if ((expr0->Type() == ET_CONSTANT && expr0->GetConstant() == 0.0) ||
        (expr1->Type() == ET_CONSTANT && expr1->GetConstant() == 0.0)) {
        return std::make_shared<Constant>(0.0);
    }
    if (expr0->Type() == ET_CONSTANT && expr0->GetConstant() == 1.0) {
        return expr1;
    }    
    if (expr0->Type() == ET_CONSTANT && expr0->GetConstant() == -1.0) {
        return -expr1;
    }    
    if (expr1->Type() == ET_CONSTANT && expr1->GetConstant() == 1.0) {
        return expr0;
    }    
    if (expr1->Type() == ET_CONSTANT && expr1->GetConstant() == -1.0) {
        return -expr0;
    }    
    if (expr0->Type() == ET_MULTIPLY && expr0->HasPartialConstant() && expr1->Type() == ET_CONSTANT) {
        auto children = expr0->Children();
        if (expr0->HasPartialConstant() == 1) {
            return (expr1->GetConstant() * children[0]->GetConstant()) * children[1];
        } else {
            return (expr1->GetConstant() * children[1]->GetConstant()) * children[0];
        }
    }
    if (expr1->Type() == ET_MULTIPLY && expr1->HasPartialConstant() && expr0->Type() == ET_CONSTANT) {
        auto children = expr1->Children();
        if (expr1->HasPartialConstant() == 1) {
            return (expr0->GetConstant() * children[0]->GetConstant()) * children[1];
        } else {
            return (expr0->GetConstant() * children[1]->GetConstant()) * children[0];
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

void EmitFunction(const std::vector<std::shared_ptr<Variable>> &inputs, 
                  const std::vector<std::shared_ptr<NamedAssignment>> &outputs,
                  const std::string &name,
                  std::ostream &os) {
    AssignmentMap assignMap;
    for (auto expr : outputs) {
        expr->Register(assignMap);
    }
    
    std::unordered_map<std::string, int> varMap;
    for (auto &var : inputs) {
        varMap[var->GetName()] = std::max(varMap[var->GetName()], var->GetIndex() + 1);
    }
    for (auto &var : outputs) {
        varMap[var->GetName()] = std::max(varMap[var->GetName()], var->GetIndex() + 1);
    }    
    
    os << "void " << name << "(";
    std::unordered_set<std::string> varSet;
    bool first = true;
    for (auto &var : inputs) {
        if (varSet.find(var->GetName()) == varSet.end()) {
            if (!first) {
                os << ", ";
            }
            first = false;
            os << "const double " << var->GetName() << "[" << varMap[var->GetName()] << "]";
            varSet.insert(var->GetName());
        }
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
    os << "\tdouble t[" << assignMap.GetAssignCount() << "];" << std::endl;
    
    for (auto expr : outputs) {
        expr->Emit(assignMap, os);        
    }
    
    os << "}" << std::endl;
}

} //namespace cdstar