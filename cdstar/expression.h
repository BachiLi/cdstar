#ifndef CDSTAR_EXPRESSION_H__
#define CDSTAR_EXPRESSION_H__

#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stack>

namespace cdstar {

class Expression;
class Boolean;
class CondExpr;

typedef std::shared_ptr<Expression> ExprPtr;
typedef std::pair<ExprPtr, ExprPtr> ExprPtrPair;

class AssignmentMap {
public:
    AssignmentMap() : m_AssignCount(0), m_BooleanCount(0), m_TabNum(1) {}
    void Register(const Expression *expr);
    void RegisterBoolean(const Expression *expr);
    int GetIndex(const Expression *expr) const;
    int GetAssignCount() const {
        return m_AssignCount;
    }
    bool IsEmitted(const Expression *expr) const {
        return m_EmittedSet.find(expr) != m_EmittedSet.end();
    }
    void SetEmitted(const Expression *expr) {
        m_EmittedSet.insert(expr);
    }
    void PushMask() {
        m_ExprMasks.push(std::unordered_set<const Expression*>());
        m_LeftCondSubtrees.push(std::unordered_set<const Expression*>());
    }
    void MaskSubtree(const Expression *expr, int rootId, bool inSubtree, bool isLeft, 
                     std::unordered_set<const Expression*> &dfsSet);
    bool IsMasked(const Expression *expr) const {
        if (m_ExprMasks.size() == 0) return false;
        return m_ExprMasks.top().find(expr) != m_ExprMasks.top().end();
    }
    void PopMask() {
        m_LeftCondSubtrees.pop();
        m_ExprMasks.pop();
    }
    void IncTab() {m_TabNum++;}
    void DecTab() {m_TabNum--;}
    std::ostream& PrintTab(std::ostream &os) const {
        for (int i = 0; i < m_TabNum; i++)
            os << "\t";
        return os;
    }
    std::vector<const CondExpr*> MergeableCondExpr(const CondExpr *condExpr) const;
private:
    bool CondPath(const CondExpr *src, const CondExpr *tgt) const {
        return m_CondPath.find({src, tgt}) != m_CondPath.end();
    }
    std::unordered_map<const Expression*, int> m_ExprMap;
    std::unordered_map<const Expression*, int> m_MaxParentId;
    std::unordered_set<const Expression*> m_EmittedSet;
    std::unordered_multimap<const Boolean*, const CondExpr*> m_CondMap;
    std::set<std::pair<const CondExpr*, const CondExpr*>> m_CondPath;
    std::stack<std::unordered_set<const Expression*>> m_ExprMasks;
    std::stack<std::unordered_set<const Expression*>> m_LeftCondSubtrees;
    int m_AssignCount;
    int m_BooleanCount;
    int m_TabNum;
};

enum ExpressionType {
    ET_VARIABLE,
    ET_CONSTANT,
    ET_INTEGER_CONSTANT,
    ET_NAMED_ASSIGNMENT,
    ET_NEGATE,
    ET_INVERSE,
    ET_SIN,
    ET_COS,
    ET_TAN,
    ET_SQRT,
    ET_ASIN,
    ET_ACOS,
    ET_LOG,
    ET_ADD,
    ET_MULTIPLY,
    ET_BOOLEAN,
    ET_CONDEXPR
};

template <typename T>
inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

class StructType;

class Argument {
public:
    Argument(const std::string &name) : m_Name(name) {}
    virtual std::string GetDeclaration() const = 0;
    virtual std::shared_ptr<Argument> SetParent(const std::string &parentName) const = 0;
    std::string GetName() const {
        return m_Name;
    }
    virtual std::shared_ptr<Argument> GetArg(const std::string &name, int index = 0) const {
        return std::shared_ptr<Argument>();
    }
    virtual std::shared_ptr<Expression> GetExpr(int index = 0) const {
        return std::shared_ptr<Expression>();
    }
    virtual std::shared_ptr<StructType> GetStructType() const {
        return std::shared_ptr<StructType>();
    }
protected:
    std::string m_Name;
};

class DoubleArgument : public Argument {
public:
    DoubleArgument(const std::string &name, int size = 1);
    std::string GetDeclaration() const;
    std::shared_ptr<Argument> SetParent(const std::string &parentName) const {
        if (parentName == "") {
            return std::make_shared<DoubleArgument>(m_Name, m_Exprs.size());
        } else {
            return std::make_shared<DoubleArgument>(parentName + "." + m_Name, m_Exprs.size());
        }
    }
    std::shared_ptr<Expression> GetExpr(int index = 0) const;
private:
    std::vector<std::shared_ptr<Expression>> m_Exprs;
};

class IntegerArgument : public Argument {
public:
    IntegerArgument(const std::string &name, int size = 1);
    std::string GetDeclaration() const;
    std::shared_ptr<Argument> SetParent(const std::string &parentName) const {
        if (parentName == "") {
            return std::make_shared<IntegerArgument>(m_Name, m_Exprs.size());
        } else {
            return std::make_shared<IntegerArgument>(parentName + "." + m_Name, m_Exprs.size());
        }
    }
    std::shared_ptr<Expression> GetExpr(int index = 0) const;
private:
    std::vector<std::shared_ptr<Expression>> m_Exprs;
};

class StructInst;

class StructType {
public:
    StructType(const std::string &name,
               const std::vector<std::shared_ptr<Argument>> &args,
               bool isUnion = false);
    std::string GetName() const {return m_Name;}
    std::shared_ptr<StructInst> GenStructInst(const std::string &name) const;
    void EmitForwardDeclaration(std::ostream &os) const;
    std::vector<std::shared_ptr<Argument>> GetArgs() const {
        return m_Args;
    }
private:
    std::string m_Name;
    std::vector<std::shared_ptr<Argument>> m_Args;
    bool m_IsUnion;
};

class StructArgument;

class StructInst {
public:
    StructInst(const std::string &name,
               const std::vector<std::shared_ptr<Argument>> &args);
    std::shared_ptr<Argument> GetArg(const std::string &argName) const {
        auto it = m_Args.find(argName);
        if (it == m_Args.end()) {
            std::cerr << "argName:" << argName << std::endl;
            throw std::runtime_error("Argument not found in StructInst");
        }
        return it->second;
    }
private:
    std::unordered_map<std::string, std::shared_ptr<Argument>> m_Args;
};

class StructArgument : public Argument {
public:
    StructArgument(const std::string &name, 
                   const StructType &structType,
                   int size = 1,
                   bool nested = false);
    std::string GetDeclaration() const;
    std::shared_ptr<Argument> SetParent(const std::string &parentName) const {
        if (parentName == "") {
            return std::make_shared<StructArgument>(
                m_Name, m_StructType, m_StructInsts.size(), true);
        } else {
            return std::make_shared<StructArgument>(
                parentName + "." + m_Name, m_StructType, m_StructInsts.size(), true);
        }
    }
    std::shared_ptr<Argument> GetArg(const std::string &argName, int index = 0) const {
        return m_StructInsts[index]->GetArg(argName);
    }
    std::shared_ptr<StructType> GetStructType() const {
        return std::make_shared<StructType>(m_StructType);
    }
private:
    StructType m_StructType;
    std::vector<std::shared_ptr<StructInst>> m_StructInsts;
    bool m_Nested;
};

class Expression {
public:
    virtual ExpressionType Type() const = 0;
    virtual void Print() const = 0;
    virtual void Register(AssignmentMap &assignMap) const;
    virtual bool UseTmpVar() const {return true;}
    virtual void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    virtual std::string GetEmitName(const AssignmentMap &assignMap) const = 0;
    virtual std::vector<std::shared_ptr<Expression>> Children() const = 0;
    virtual std::vector<std::shared_ptr<Expression>> Dervs() const = 0;
    virtual int HasPartialConstant() const {return 0;}
    virtual double GetConstant() const {return 0.0;}    
    virtual bool Equal(const std::shared_ptr<Expression> expr) const = 0;
    size_t GetHash() const {return m_Hash;}
protected:
    virtual void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const = 0;
    size_t m_Hash;
};

class Variable : public Expression {
public:
    Variable(const std::string &name, const int index) : 
            m_Name(name), m_Index(index) {
        m_Hash = ComputeHash();
    }
    bool UseTmpVar() const {return false;}
    ExpressionType Type() const {return ET_VARIABLE;}
    void Print() const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::string GetName() const {return m_Name;}
    int GetIndex() const {return m_Index;}
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<Variable> _expr = std::dynamic_pointer_cast<Variable>(expr);
        return m_Name == _expr->m_Name && m_Index == _expr->m_Index;
    }
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(m_Index);
        hash_combine(hash, std::hash<std::string>()(m_Name));
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
private:
    std::string m_Name;  
    int m_Index;
};

class Constant : public Expression {
public:
    Constant(const double value) : m_Value(value) {
        m_Hash = ComputeHash();
    }
    bool UseTmpVar() const {return false;}
    ExpressionType Type() const {return ET_CONSTANT;}
    void Print() const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;   
    double GetConstant() const {return m_Value;}
    std::shared_ptr<Expression> ReplaceConstant(double v) {
        return std::make_shared<Constant>(v);
    }
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<Constant> _expr = std::dynamic_pointer_cast<Constant>(expr);
        return m_Value == _expr->m_Value;
    }
    size_t ComputeHash() const {
        std::size_t hash = std::hash<double>()(m_Value);
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;    
private:
    double m_Value;
};

class IntegerConstant : public Expression {
public:
    IntegerConstant(const int value) : m_Value(value) {
        m_Hash = ComputeHash();
    }
    bool UseTmpVar() const {return false;}
    ExpressionType Type() const {return ET_INTEGER_CONSTANT;}
    void Print() const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<IntegerConstant> _expr = std::dynamic_pointer_cast<IntegerConstant>(expr);
        return m_Value == _expr->m_Value;
    }
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(m_Value);
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;    
private:
    int m_Value;
};

class NamedAssignment : public Expression {
public:
    NamedAssignment(const std::string &name, const int index,
                    const std::shared_ptr<Expression> &expr) :
            m_Name(name), m_Index(index), m_Expr(expr) {
        m_Hash = ComputeHash();
    }
    bool UseTmpVar() const {return false;}
    ExpressionType Type() const {return ET_NAMED_ASSIGNMENT;}
    void Print() const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::string GetName() const {return m_Name;}
    int GetIndex() const {return m_Index;}    
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<NamedAssignment> _expr = std::dynamic_pointer_cast<NamedAssignment>(expr);
        return m_Name == _expr->m_Name && m_Expr->Equal(_expr->m_Expr);
    }    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(m_Index);
        hash_combine(hash, std::hash<std::string>()(m_Name));
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;    
private:    
    std::string m_Name;
    int m_Index;
    std::shared_ptr<Expression> m_Expr;
};

class UnaryAssignment : public Expression {
public:
    UnaryAssignment(const std::shared_ptr<Expression> &expr) :
        m_Expr(expr) {}
    virtual std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;
    virtual bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<UnaryAssignment> _expr = std::dynamic_pointer_cast<UnaryAssignment>(expr);
        return m_Expr->Equal(_expr->m_Expr);
    }
protected:
    std::shared_ptr<Expression> m_Expr;
};

class Negate : public UnaryAssignment {
public:
    Negate(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_NEGATE;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_NEGATE);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Inverse : public UnaryAssignment, public std::enable_shared_from_this<Inverse> {
public:
    Inverse(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_INVERSE;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_INVERSE);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;    
};

class Sin : public UnaryAssignment {
public:
    Sin(const std::shared_ptr<Expression> &expr) :
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_SIN;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_SIN);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Cos : public UnaryAssignment {
public:
    Cos(const std::shared_ptr<Expression> &expr) :
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_COS;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_COS);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Tan : public UnaryAssignment {
public:
    Tan(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_TAN;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_TAN);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Sqrt : public UnaryAssignment, public std::enable_shared_from_this<Sqrt> {
public:
    Sqrt(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_SQRT;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_SQRT);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class ASin : public UnaryAssignment {
public:
    ASin(const std::shared_ptr<Expression> &expr) :
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_ASIN;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_ASIN);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class ACos : public UnaryAssignment {
public:
    ACos(const std::shared_ptr<Expression> &expr) :
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_ACOS;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_ACOS);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Log : public UnaryAssignment {
public:
    Log(const std::shared_ptr<Expression> &expr) :
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_LOG;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_LOG);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class BinaryAssignment : public Expression {
public:
    BinaryAssignment(const std::shared_ptr<Expression> expr0,
                     const std::shared_ptr<Expression> expr1) : 
        m_Expr0(expr0), m_Expr1(expr1) {
    }
    virtual std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;  
    int HasPartialConstant() const {
        if (m_Expr0->Type() == ET_CONSTANT) {
            return 1;
        } else if (m_Expr1->Type() == ET_CONSTANT) {
            return 2;
        } else {
            return 0;
        }
    }
protected:
    std::shared_ptr<Expression> m_Expr0, m_Expr1;
};

class Add : public BinaryAssignment {
public:    
    Add(const std::shared_ptr<Expression> expr0,
        const std::shared_ptr<Expression> expr1) : 
            BinaryAssignment(expr0, expr1) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_ADD;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<Add> _expr = std::dynamic_pointer_cast<Add>(expr);
        return m_Expr0->Equal(_expr->m_Expr0) && m_Expr1->Equal(_expr->m_Expr1);
    }    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_ADD);
        hash_combine(hash, m_Expr0->GetHash());
        hash_combine(hash, m_Expr1->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Multiply : public BinaryAssignment {
public:    
    Multiply(const std::shared_ptr<Expression> expr0,
             const std::shared_ptr<Expression> expr1) :
                BinaryAssignment(expr0, expr1) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_MULTIPLY;}
    void Print() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;   
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<Multiply> _expr = std::dynamic_pointer_cast<Multiply>(expr);
        return m_Expr0->Equal(_expr->m_Expr0) && m_Expr1->Equal(_expr->m_Expr1);
    }        
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_MULTIPLY);
        hash_combine(hash, m_Expr0->GetHash());
        hash_combine(hash, m_Expr1->GetHash());
        return hash;
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
};

class Boolean : public Expression {
public:
    enum Op {
        GREATER,
        GREATER_OR_EQUAL,
        EQUAL,
        NOT_EQUAL,
        LESS_OR_EQUAL,
        LESS
    };
    Boolean(const Op op,
            const std::shared_ptr<Expression> expr0,
            const std::shared_ptr<Expression> expr1) :
            m_Op(op), m_Expr0(expr0), m_Expr1(expr1) {
        m_Hash = ComputeHash();
    }    
    ExpressionType Type() const {return ET_BOOLEAN;}
    void Print() const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<Boolean> _expr = std::dynamic_pointer_cast<Boolean>(expr);        
        return m_Op == _expr->m_Op && m_Expr0->Equal(_expr->m_Expr0) && m_Expr1->Equal(_expr->m_Expr1);
    }    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(Type());
        hash_combine(hash, std::hash<int>()(m_Op));
        hash_combine(hash, m_Expr0->GetHash());
        hash_combine(hash, m_Expr1->GetHash());
        return hash;
    }
    std::shared_ptr<Boolean> Swap() const {
        return std::make_shared<Boolean>(ReverseOp(), m_Expr1, m_Expr0);
    }
protected:    
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
private:
    std::string OpToString() const;
    Op ReverseOp() const;
    Op m_Op;
    std::shared_ptr<Expression> m_Expr0, m_Expr1;
};

class CondExpr : public Expression {
public:
    CondExpr(const std::shared_ptr<Boolean> cond,
             const std::shared_ptr<Expression> trueExpr,
             const std::shared_ptr<Expression> falseExpr) :
             m_Cond(cond), m_TrueExpr(trueExpr), m_FalseExpr(falseExpr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_CONDEXPR;}
    void Print() const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::string GetEmitName(const AssignmentMap &assignMap) const;
    std::vector<std::shared_ptr<Expression>> Children() const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    bool Equal(const std::shared_ptr<Expression> expr) const {
        if (expr->Type() != Type() || GetHash() != expr->GetHash()) return false;
        std::shared_ptr<CondExpr> _expr = std::dynamic_pointer_cast<CondExpr>(expr);
        return m_Cond->Equal(_expr->m_Cond) &&
               m_TrueExpr->Equal(_expr->m_TrueExpr) &&
               m_FalseExpr->Equal(_expr->m_FalseExpr);
    }
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(Type());
        hash_combine(hash, m_Cond->GetHash());
        hash_combine(hash, m_TrueExpr->GetHash());
        hash_combine(hash, m_FalseExpr->GetHash());
        return hash;
    }
    std::shared_ptr<Boolean> GetCond() const {
        return m_Cond;
    }
    std::shared_ptr<CondExpr> Swap() const {
        return std::make_shared<CondExpr>(m_Cond->Swap(), m_FalseExpr, m_TrueExpr);
    }
protected:
    void EmitSelf(AssignmentMap &assignMap, std::ostream &os) const;
private:
    std::shared_ptr<Boolean> m_Cond;
    std::shared_ptr<Expression> m_TrueExpr, m_FalseExpr;
};

void PrintExpressionCache();
void ClearExpressionCache();
std::shared_ptr<Expression> CacheExpression(const std::shared_ptr<Expression> &expr);
std::shared_ptr<Expression> CacheExpression(const std::shared_ptr<Expression> &expr0,
                                            const std::shared_ptr<Expression> &expr1);

std::shared_ptr<Expression> operator+(const std::shared_ptr<Expression> expr0,
                                      const std::shared_ptr<Expression> expr1);
std::shared_ptr<Expression> operator+(const double v0,
                                      const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Expression> operator+(const std::shared_ptr<Expression> expr0,
                                             const double v1) {
    return v1 + expr0;
}
inline std::shared_ptr<Expression>& operator+=(std::shared_ptr<Expression> &expr0, 
                                               const std::shared_ptr<Expression> &expr1) {
    expr0 = expr0 + expr1;
    return expr0;
}
inline std::shared_ptr<Expression>& operator+=(std::shared_ptr<Expression> &expr0, 
                                               const double v1) {
    expr0 = expr0 + v1;
    return expr0;
}

std::shared_ptr<Expression> operator-(const std::shared_ptr<Expression> expr);
std::shared_ptr<Expression> operator-(const std::shared_ptr<Expression> expr0, 
                                      const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Expression> operator-(const double v0,
                                             const std::shared_ptr<Expression> expr1) {
    return v0 + (-expr1);
}
inline std::shared_ptr<Expression> operator-(const std::shared_ptr<Expression> expr0,
                                             const double v1) {
    return expr0 + (-v1);
}
inline std::shared_ptr<Expression>& operator-=(std::shared_ptr<Expression> &expr0, 
                                               const std::shared_ptr<Expression> &expr1) {
    expr0 = expr0 - expr1;
    return expr0;
}
inline std::shared_ptr<Expression>& operator-=(std::shared_ptr<Expression> &expr0, 
                                               const double v1) {
    expr0 = expr0 - v1;
    return expr0;
}

std::shared_ptr<Expression> operator*(const std::shared_ptr<Expression> expr0,
                                      const std::shared_ptr<Expression> expr1);
std::shared_ptr<Expression> operator*(const double v0,
                                      const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Expression> operator*(const std::shared_ptr<Expression> expr0,
                                             const double v1) {
    return v1 * expr0;
}
inline std::shared_ptr<Expression>& operator*=(std::shared_ptr<Expression> &expr0, 
                                               const std::shared_ptr<Expression> &expr1) {
    expr0 = expr0 * expr1;
    return expr0;
}
inline std::shared_ptr<Expression>& operator*=(std::shared_ptr<Expression> &expr0, 
                                               const double v1) {
    expr0 = expr0 * v1;
    return expr0;
}

std::shared_ptr<Expression> Inv(const std::shared_ptr<Expression> expr);

std::shared_ptr<Expression> operator/(const std::shared_ptr<Expression> expr0,
                                      const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Expression> operator/(const double v0,
                                             const std::shared_ptr<Expression> expr1) {
    return v0 * Inv(expr1);
}
inline std::shared_ptr<Expression> operator/(const std::shared_ptr<Expression> expr0,
                                             const double v1) {
    return expr0 * (1.0 / v1);
}
inline std::shared_ptr<Expression>& operator/=(std::shared_ptr<Expression> &expr0, 
                                               const std::shared_ptr<Expression> &expr1) {
    expr0 = expr0 / expr1;
    return expr0;
}
inline std::shared_ptr<Expression>& operator/=(std::shared_ptr<Expression> &expr0, 
                                               const double v1) {
    expr0 = expr0 / v1;
    return expr0;
}

inline std::shared_ptr<Expression> sin(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Sin>(expr));
}

inline std::shared_ptr<Expression> cos(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Cos>(expr));
}

inline std::shared_ptr<Expression> tan(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Tan>(expr));
}

inline std::shared_ptr<Expression> asin(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<ASin>(expr));
}

inline std::shared_ptr<Expression> acos(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<ACos>(expr));
}

inline std::shared_ptr<Expression> log(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Log>(expr));
}

inline std::shared_ptr<Expression> sqrt(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Sqrt>(expr));
}

// >
std::shared_ptr<Boolean> Gt(const std::shared_ptr<Expression> expr0,
                            const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Gt(const double v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Gt(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Gt(const std::shared_ptr<Expression> expr0,
                                   const double v1) {
    return Gt(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Gt(const int v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Gt(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Gt(const std::shared_ptr<Expression> expr0,
                                   const int v1) {
    return Gt(expr0, std::make_shared<IntegerConstant>(v1));
}

// >=
std::shared_ptr<Boolean> Gte(const std::shared_ptr<Expression> expr0,
                             const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Gte(const double v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Gte(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Gte(const std::shared_ptr<Expression> expr0,
                                    const double v1) {
    return Gte(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Gte(const int v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Gte(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Gte(const std::shared_ptr<Expression> expr0,
                                    const int v1) {
    return Gte(expr0, std::make_shared<IntegerConstant>(v1));
}
 
// ==                                   
std::shared_ptr<Boolean> Eq(const std::shared_ptr<Expression> expr0,
                            const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Eq(const double v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Eq(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Eq(const std::shared_ptr<Expression> expr0,
                                   const double v1) {
    return Eq(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Eq(const int v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Eq(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Eq(const std::shared_ptr<Expression> expr0,
                                   const int v1) {
    return Eq(expr0, std::make_shared<IntegerConstant>(v1));
}

// !=
std::shared_ptr<Boolean> Neq(const std::shared_ptr<Expression> expr0,
                                    const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Neq(const double v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Neq(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Neq(const std::shared_ptr<Expression> expr0,
                                    const double v1) {
    return Neq(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Neq(const int v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Neq(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Neq(const std::shared_ptr<Expression> expr0,
                                    const int v1) {
    return Neq(expr0, std::make_shared<IntegerConstant>(v1));
}

// <=
std::shared_ptr<Boolean> Lte(const std::shared_ptr<Expression> expr0,
                             const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Lte(const double v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Lte(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Lte(const std::shared_ptr<Expression> expr0,
                                    const double v1) {
    return Lte(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Lte(const int v0,
                                    const std::shared_ptr<Expression> expr1) {
    return Lte(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Lte(const std::shared_ptr<Expression> expr0,
                                    const int v1) {
    return Lte(expr0, std::make_shared<IntegerConstant>(v1));
}

// <
std::shared_ptr<Boolean> Lt(const std::shared_ptr<Expression> expr0,
                                   const std::shared_ptr<Expression> expr1);
inline std::shared_ptr<Boolean> Lt(const double v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Lt(std::make_shared<Constant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Lt(const std::shared_ptr<Expression> expr0,
                                   const double v1) {
    return Lt(expr0, std::make_shared<Constant>(v1));
}
inline std::shared_ptr<Boolean> Lt(const int v0,
                                   const std::shared_ptr<Expression> expr1) {
    return Lt(std::make_shared<IntegerConstant>(v0), expr1);
}
inline std::shared_ptr<Boolean> Lt(const std::shared_ptr<Expression> expr0,
                                   const int v1) {
    return Lt(expr0, std::make_shared<IntegerConstant>(v1));
}

std::shared_ptr<Expression> IfElse(const std::shared_ptr<Boolean> cond,
                                   const std::shared_ptr<Expression> trueExpr,
                                   const std::shared_ptr<Expression> falseExpr);
inline std::shared_ptr<Expression> IfElse(const std::shared_ptr<Boolean> cond,
                                          const std::shared_ptr<Expression> trueExpr,
                                          const double falseExpr) {
    return IfElse(cond, trueExpr, std::make_shared<Constant>(falseExpr));
}
inline std::shared_ptr<Expression> IfElse(const std::shared_ptr<Boolean> cond,
                                          const double trueExpr,
                                          const std::shared_ptr<Expression> falseExpr) {
    return IfElse(cond, std::make_shared<Constant>(trueExpr), falseExpr);
}

std::vector<std::shared_ptr<Expression>> Derivatives(const std::vector<ExprPtrPair> &dervExprs);

void EmitFunction(const std::vector<std::shared_ptr<Argument>> &input,
                  const std::vector<std::shared_ptr<NamedAssignment>> &output,
                  const std::string &name,
                  std::ostream &os);

} //namespace cdstar

#endif //CDSTAR_EXPRESSION_H__