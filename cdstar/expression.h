#ifndef CDSTAR_EXPRESSION_H__
#define CDSTAR_EXPRESSION_H__

#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace cdstar {

class Expression;

typedef std::shared_ptr<Expression> ExprPtr;
typedef std::pair<ExprPtr, ExprPtr> ExprPtrPair;

class AssignmentMap {
public:
    AssignmentMap() : m_AssignCount(0), m_BooleanCount(0) {}
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
private:
    std::unordered_map<const Expression*, int> m_ExprMap;
    std::unordered_set<const Expression*> m_EmittedSet;
    int m_AssignCount;
    int m_BooleanCount;
};

enum ExpressionType {
    ET_VARIABLE,
    ET_CONSTANT,
    ET_NAMED_ASSIGNMENT,    
    ET_NEGATE,
    ET_INVERSE,
    ET_SIN,
    ET_COS,
    ET_TAN,
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

class Expression {
public:
    virtual ExpressionType Type() const = 0;
    virtual void Print() const = 0;
    virtual void Register(AssignmentMap &assignMap) const = 0;
    virtual void Emit(AssignmentMap &assignMap, std::ostream &os) const = 0;
    virtual std::string GetEmitName(const AssignmentMap &assignMap) const = 0;
    virtual std::vector<std::shared_ptr<Expression>> Children() const = 0;
    virtual std::vector<std::shared_ptr<Expression>> Dervs() const = 0;
    virtual int HasPartialConstant() const {return 0;}
    virtual double GetConstant() const {return 0.0;}    
    virtual bool Equal(const std::shared_ptr<Expression> expr) const = 0;
    size_t GetHash() const {return m_Hash;}
protected:
    size_t m_Hash;
};

class Variable : public Expression {
public:
    Variable(const std::string &name, const int index) : 
            m_Name(name), m_Index(index) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_VARIABLE;}
    void Print() const;
    void Register(AssignmentMap &assignMap) const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
private:
    std::string m_Name;  
    int m_Index;
};

class Constant : public Expression {
public:
    Constant(const double value) : m_Value(value) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_CONSTANT;}
    void Print() const;
    void Register(AssignmentMap &assignMap) const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
private:
    double m_Value;
};

class NamedAssignment : public Expression {
public:
    NamedAssignment(const std::string &name, const int index,
                    const std::shared_ptr<Expression> &expr) :
            m_Name(name), m_Index(index), m_Expr(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_NAMED_ASSIGNMENT;}
    void Print() const;
    void Register(AssignmentMap &assignMap) const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
private:    
    std::string m_Name;
    int m_Index;
    std::shared_ptr<Expression> m_Expr;
};

class UnaryAssignment : public Expression {
public:
    UnaryAssignment(const std::shared_ptr<Expression> &expr) :
        m_Expr(expr) {}
    virtual void Register(AssignmentMap &assignMap) const;
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
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_NEGATE);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }    
};

class Inverse : public UnaryAssignment {
public:
    Inverse(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_INVERSE;}
    void Print() const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_INVERSE);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }    
};

class Sin : public UnaryAssignment {
public:
    Sin(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_SIN;}
    void Print() const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_SIN);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }    
};

class Cos : public UnaryAssignment {
public:
    Cos(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_COS;}
    void Print() const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_COS);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }    
};

class Tan : public UnaryAssignment {
public:
    Tan(const std::shared_ptr<Expression> &expr) : 
            UnaryAssignment(expr) {
        m_Hash = ComputeHash();
    }
    ExpressionType Type() const {return ET_TAN;}
    void Print() const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
    std::vector<std::shared_ptr<Expression>> Dervs() const;    
    size_t ComputeHash() const {
        std::size_t hash = std::hash<int>()(ET_TAN);
        hash_combine(hash, m_Expr->GetHash());
        return hash;
    }
};

class BinaryAssignment : public Expression {
public:
    BinaryAssignment(const std::shared_ptr<Expression> expr0,
                     const std::shared_ptr<Expression> expr1) : 
        m_Expr0(expr0), m_Expr1(expr1) {
    }
    virtual void Register(AssignmentMap &assignMap) const;
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
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
    void Register(AssignmentMap &assignMap) const;
    void Emit(AssignmentMap &assignMap, std::ostream &os) const;
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
    void Register(AssignmentMap &assignMap) const;
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
    std::shared_ptr<CondExpr> Swap() const {
        return std::make_shared<CondExpr>(m_Cond->Swap(), m_FalseExpr, m_TrueExpr);
    }    
private:
    std::shared_ptr<Boolean> m_Cond;
    std::shared_ptr<Expression> m_TrueExpr, m_FalseExpr;
};

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

std::shared_ptr<Expression> operator*(const std::shared_ptr<Expression> expr0, 
                                      const std::shared_ptr<Expression> expr1);

std::shared_ptr<Expression> operator*(const double v0, 
                                      const std::shared_ptr<Expression> expr1);

inline std::shared_ptr<Expression> operator*(const std::shared_ptr<Expression> expr0, 
                                             const double v1) {
    return v1 * expr0;
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

inline std::shared_ptr<Expression> sin(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Sin>(expr));
}

inline std::shared_ptr<Expression> cos(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Cos>(expr));
}

inline std::shared_ptr<Expression> tan(const std::shared_ptr<Expression> expr) {
    return CacheExpression(std::make_shared<Tan>(expr));
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

std::shared_ptr<Expression> IfElse(const std::shared_ptr<Boolean> cond,
                                   const std::shared_ptr<Expression> trueExpr,
                                   const std::shared_ptr<Expression> falseExpr);

std::vector<std::shared_ptr<Expression>> Derivatives(const std::vector<ExprPtrPair> &dervExprs);

void EmitFunction(const std::vector<std::shared_ptr<Variable>> &input,
                  const std::vector<std::shared_ptr<NamedAssignment>> &output,
                  const std::string &name,
                  std::ostream &os);

} //namespace cdstar

#endif //CDSTAR_EXPRESSION_H__