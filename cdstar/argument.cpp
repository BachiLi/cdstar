#include "argument.h"
#include "expression.h"

namespace cdstar {

DoubleArgument::DoubleArgument(const std::string &name, int size) :
        Argument(name) {
    for (int i = 0; i < size; i++) {
        m_Exprs.push_back(std::make_shared<Variable>(name, size == 1 ? -1 : i));
    }
}

std::string DoubleArgument::GetDeclaration() const {
    if (m_Exprs.size() == 1) {
        return "double " + m_Name;
    } else {
        return "double " + m_Name + "[" + std::to_string(m_Exprs.size()) + "]";
    }
}

std::shared_ptr<Expression> DoubleArgument::GetExpr(int index) const {
    return m_Exprs[index];
}

IntegerArgument::IntegerArgument(const std::string &name, int size) :
        Argument(name) {
    for (int i = 0; i < size; i++) {
        m_Exprs.push_back(std::make_shared<Variable>(name, size == 1 ? -1 : i));
    }
}

std::string IntegerArgument::GetDeclaration() const {
    if (m_Exprs.size() == 1) {
        return "int " + m_Name;
    } else {
        return "int " + m_Name + "[" + std::to_string(m_Exprs.size()) + "]";
    }
}

std::shared_ptr<Expression> IntegerArgument::GetExpr(int index) const {
    return m_Exprs[index];
}

StructType::StructType(const std::string &name, 
                       const std::vector<std::shared_ptr<Argument>> &args,
                       bool isUnion) :
        m_Name(name), m_IsUnion(isUnion) {
    for (auto &arg : args) {
        m_Args.push_back(arg->SetParent(""));
    }
}

std::shared_ptr<StructInst> StructType::GenStructInst(const std::string &name) const {
    return std::make_shared<StructInst>(name, m_Args);
}

void StructType::EmitForwardDeclaration(std::ostream &os) const {
    if (!m_IsUnion) {
        os << "typedef struct {" << std::endl;
    } else {
        os << "typedef union {" << std::endl;
    }
    for (auto &arg : m_Args) {
        os << "\t" << arg->GetDeclaration() << ";" << std::endl;
    }
    os << "} " << m_Name << ";" << std::endl;
}

StructInst::StructInst(const std::string &name,
                       const std::vector<std::shared_ptr<Argument>> &args) {
    for (auto &arg : args) {
        m_Args[arg->GetName()] = arg->SetParent(name);
    }
}

StructArgument::StructArgument(const std::string &name, 
                               const StructType &structType,
                               int size,
                               bool nested) :
        Argument(name), m_StructType(structType), m_Nested(nested) {
    for (int i = 0; i < size; i++) {
        if (m_Nested && size == 1) {
            m_StructInsts.push_back(m_StructType.GenStructInst(m_Name));
        } else {
            m_StructInsts.push_back(m_StructType.GenStructInst(
                m_Name + "[" + std::to_string(i) + "]"));
        }
    }
}

std::string StructArgument::GetDeclaration() const {
    if (m_Nested && m_StructInsts.size() == 1) {
        return m_StructType.GetName() + " " + m_Name;
    } else {
        return m_StructType.GetName() + " " + m_Name + "[" + std::to_string(m_StructInsts.size()) + "]";
    }
}

}