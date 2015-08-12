#ifndef CDSTAR_ARGUMENT_H__
#define CDSTAR_ARGUMENT_H__

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

namespace cdstar {

class StructType;
class Expression;

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

}

#endif //CDSTAR_ARGUMENT_H__