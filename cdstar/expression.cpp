#include "expression.h"
#include "derivativegraph.h"

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
        return expr0->Equal(expr1);
    }
};

std::unordered_set<std::shared_ptr<Expression>, ExprHasher, ExprComparator> g_ExprCache;

void AssignmentMap::Register(const Expression *expr) {
    if (m_ExprMap.find(expr) == m_ExprMap.end()) {
        m_ExprMap[expr] = m_AssignCount++;
        m_MaxParentId[expr] = -1;
        
        if (expr->Type() == ET_CONDEXPR) {
            const CondExpr *_expr = dynamic_cast<const CondExpr*>(expr);
            m_CondMap.insert({_expr->GetCond().get(), _expr});
            std::stack<const Expression*> dfsStack;
            std::unordered_set<const Expression*> dfsSet;
            dfsStack.push(expr);
            while (!dfsStack.empty()) {
                const Expression *curExpr = dfsStack.top();
                dfsStack.pop();
                if (dfsSet.find(curExpr) == dfsSet.end()) {
                    continue;
                }
                dfsSet.insert(curExpr);
                if (curExpr != expr && curExpr->Type() == ET_CONDEXPR) {
                    const CondExpr *_curExpr = dynamic_cast<const CondExpr*>(curExpr);
                    m_CondPath.insert({_expr, _curExpr});
                }
                for (auto child : curExpr->Children()) {
                    dfsStack.push(child.get());
                }
            }
        }
    }
    int id = m_ExprMap[expr];
    for (auto child : expr->Children()) {
        m_MaxParentId[child.get()] = std::max(m_MaxParentId[child.get()], id);
    }
}

void AssignmentMap::RegisterBoolean(const Expression *expr) {
    if (m_ExprMap.find(expr) == m_ExprMap.end()) {
        m_ExprMap[expr] = m_BooleanCount++;
    }
}

int AssignmentMap::GetIndex(const Expression *expr) const {
    auto it = m_ExprMap.find(expr);
    if (it != m_ExprMap.end()) {
        return it->second;
    }
    return -1;
}

void AssignmentMap::MaskSubtree(const Expression *expr, int rootId, 
                                bool inSubtree, bool isLeft) {
    if (IsEmitted(expr)) {
        return;        
    }
    bool out = rootId < 0 || m_MaxParentId[expr] > rootId;
    if ((inSubtree && !out) || (!inSubtree && out)) {  
        m_ExprMasks.top().insert(expr);
    }
    
    if (isLeft) {
        m_LeftCondSubtrees.top().insert(expr);
    } else {
        if (m_LeftCondSubtrees.top().find(expr) != m_LeftCondSubtrees.top().end()) {
            if (inSubtree) {
                // Exception: the intersection of the left and right trees should
                //            be unmasked
                m_ExprMasks.top().erase(expr);
            } else {
                // Exception: the intersection of the left and right trees should
                //            also be masked                
                m_ExprMasks.top().insert(expr);
            }
        }
    }
    
    if (out) {
        rootId = -1;
    }
    for (auto child : expr->Children()) {
        MaskSubtree(child.get(), rootId, inSubtree, isLeft);
    }
}

std::vector<const CondExpr*> AssignmentMap::MergeableCondExpr(const CondExpr *condExpr) const {
    std::vector<const CondExpr*> ret;
    auto ranges = m_CondMap.equal_range(condExpr->GetCond().get());
    for (auto it = ranges.first; it != ranges.second; it++) {
        const CondExpr *srcExpr = it->second;
        if (m_CondPath.find({srcExpr, condExpr}) == m_CondPath.end()) {
            ret.push_back(srcExpr);
        }
    }
    return ret;
}

DoubleArgument::DoubleArgument(const std::string &name) :
    m_Name(name), m_Expr(std::make_shared<Variable>(name, -1)) {}

std::string DoubleArgument::GetDeclaration() const {
    return "double " + m_Name;
}

std::shared_ptr<Expression> DoubleArgument::GetExpr() const {
    return m_Expr;
}

DoubleArrayArgument::DoubleArrayArgument(const std::string &name, int size) :
        m_Name(name) {
    for (int i = 0; i < size; i++) {
        m_Exprs.push_back(std::make_shared<Variable>(name, i));
    }
}

std::string DoubleArrayArgument::GetDeclaration() const {
    return "double " + m_Name + "[" + std::to_string(m_Exprs.size()) + "]";
}

std::shared_ptr<Expression> DoubleArrayArgument::GetExpr(int index) const {
    return m_Exprs[index];
}

void Expression::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    auto children = Children();
    if (children.size() == 0) {
        return;
    }
    for (auto child : children) {
        child->Emit(assignMap, os);
    }
    if (!assignMap.IsMasked(this)) {
        EmitSelf(assignMap, os);
        assignMap.SetEmitted(this);
    }
}

void Expression::Register(AssignmentMap &assignMap) const {
    auto children = Children();
    if (children.size() == 0) {
        return;
    }
    for (auto child : children) {
        child->Register(assignMap);
    }
    if (UseTmpVar()) {
        assignMap.Register(this);
    }
}

void Variable::Print() const {
    if (m_Index >= 0) {
        std::cerr << m_Name << "[" << m_Index << "]";
    } else {
        std::cerr << m_Name << std::endl;
    }
}

void Variable::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
}

std::string Variable::GetEmitName(const AssignmentMap &assignMap) const {
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

void Constant::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
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

void NamedAssignment::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = " << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
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

void Negate::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = -" << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Negate::Dervs() const {
    return {std::make_shared<Constant>(-1.0)};
}

void Inverse::Print() const {
    std::cerr << "(1.0 / ";
    m_Expr->Print();
    std::cerr << ")";
}

void Inverse::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = 1.0 / " << m_Expr->GetEmitName(assignMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Inverse::Dervs() const {
    // d(1/x)/dx = -1/x^2
    std::shared_ptr<Inverse> thisPtr = std::const_pointer_cast<Inverse>(shared_from_this());
    return {-(thisPtr * thisPtr)};
}

void Sin::Print() const {
    std::cerr << "sin(";
    m_Expr->Print();
    std::cerr << ")";
}

void Sin::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = sin(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Sin::Dervs() const {    
    return {cos(m_Expr)};
}

void Cos::Print() const {
    std::cerr << "cos(";
    m_Expr->Print();
    std::cerr << ")";
}

void Cos::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = cos(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Cos::Dervs() const {    
    return {-sin(m_Expr)};
}

void Tan::Print() const {
    std::cerr << "tan(";
    m_Expr->Print();
    std::cerr << ")";
}

void Tan::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = tan(" << m_Expr->GetEmitName(assignMap) << ");" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Tan::Dervs() const {
    auto c = cos(m_Expr);
    return {Inv(c * c)};
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

void Add::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = " << 
        m_Expr0->GetEmitName(assignMap) << " + " <<
        m_Expr1->GetEmitName(assignMap) << ";" << std::endl;
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

void Multiply::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << GetEmitName(assignMap) << " = " << 
        m_Expr0->GetEmitName(assignMap) << " * " << 
        m_Expr1->GetEmitName(assignMap) << ";" << std::endl;
}

std::vector<std::shared_ptr<Expression>> Multiply::Dervs() const {
    return {m_Expr1, m_Expr0};
}    

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

void Boolean::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
    assignMap.PrintTab(os) << "int " << GetEmitName(assignMap) << " = " <<
        m_Expr0->GetEmitName(assignMap) << OpToString() << m_Expr1->GetEmitName(assignMap) << ";" << std::endl;
}

std::string Boolean::GetEmitName(const AssignmentMap &assignMap) const {
    return std::string("b") + std::to_string(assignMap.GetIndex(this));
}

std::vector<std::shared_ptr<Expression>> Boolean::Children() const {
    return {m_Expr0, m_Expr1};
}

std::vector<std::shared_ptr<Expression>> Boolean::Dervs() const {
    auto zero = std::make_shared<Constant>(0.0);
    return {zero, zero};
}

void CondExpr::Print() const {
    std::cerr << "(";
    m_Cond->Print();
    std::cerr << " ? ";
    m_TrueExpr->Print();
    std::cerr << " : ";
    m_FalseExpr->Print();
    std::cerr << ")";
}

void CondExpr::Emit(AssignmentMap &assignMap, std::ostream &os) const {
    if (assignMap.IsEmitted(this)) {
        return;
    }
    m_Cond->Emit(assignMap, os);
    if (assignMap.IsMasked(this)) {
        m_TrueExpr->Emit(assignMap, os);
        m_FalseExpr->Emit(assignMap, os);
        return;
    }
    // Gather condition expressions to merge
    auto condExprs = assignMap.MergeableCondExpr(this);
    std::sort(condExprs.begin(), condExprs.end(), 
        [&](const CondExpr *e0, const CondExpr *e1) {
            return assignMap.GetIndex(e0) < assignMap.GetIndex(e1);
        });
    assignMap.PushMask();
    int id = -1;
    for (auto expr : condExprs) {
        id = std::max(id, assignMap.GetIndex(expr));
    }
    for (auto expr : condExprs) {
        assignMap.MaskSubtree(expr->m_TrueExpr.get(), id, true, true);        
    }
    for (auto expr : condExprs) {
        assignMap.MaskSubtree(expr->m_FalseExpr.get(), id, true, false);
    }
    for (auto expr : condExprs) {
        expr->m_TrueExpr->Emit(assignMap, os);
        expr->m_FalseExpr->Emit(assignMap, os);
    }
    assignMap.PopMask();
    assignMap.PrintTab(os) << "if (" << m_Cond->GetEmitName(assignMap) << ") {" << std::endl;
    assignMap.IncTab();
    assignMap.PushMask();
    for (auto expr : condExprs) {
        assignMap.MaskSubtree(expr->m_TrueExpr.get(), id, false, true);
    }
    for (auto expr : condExprs) {
        assignMap.MaskSubtree(expr->m_FalseExpr.get(), id, false, false);
    }
    for (auto expr : condExprs) {
        expr->m_TrueExpr->Emit(assignMap, os);
    }
    for (auto expr : condExprs) {
        assignMap.PrintTab(os) << expr->GetEmitName(assignMap) << " = " <<
            expr->m_TrueExpr->GetEmitName(assignMap) << ";" << std::endl;
    }
    assignMap.DecTab();
    assignMap.PrintTab(os) << "} else {" << std::endl;
    assignMap.IncTab();    
    for (auto expr : condExprs) {
        expr->m_FalseExpr->Emit(assignMap, os);
    }
    assignMap.PopMask();
    for (auto expr : condExprs) {
        assignMap.PrintTab(os) << expr->GetEmitName(assignMap) << " = " <<
            expr->m_FalseExpr->GetEmitName(assignMap) << ";" << std::endl;
    }
    assignMap.DecTab();
    assignMap.PrintTab(os) << "}" << std::endl;
    for (auto expr : condExprs) {
        assignMap.SetEmitted(expr);
    }
}

void CondExpr::EmitSelf(AssignmentMap &assignMap, std::ostream &os) const {
}

std::string CondExpr::GetEmitName(const AssignmentMap &assignMap) const {
    return std::string("t[") + std::to_string(assignMap.GetIndex(this)) + std::string("]");
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
    AssignmentMap assignMap;
    for (auto expr : outputs) {
        expr->Register(assignMap);
    }

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
    os << "\tdouble t[" << assignMap.GetAssignCount() << "];" << std::endl;
    
    for (auto expr : outputs) {
        expr->Emit(assignMap, os);
    }
    
    os << "}" << std::endl;
}

} //namespace cdstar