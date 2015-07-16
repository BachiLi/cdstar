#ifndef CDSTAR_LIBRARY_H__
#define CDSTAR_LIBRARY_H__

#include <memory>
#include <fstream>
#include <vector>
#include <unordered_set>

#ifdef _WIN32
    #include <windows.h>
    typedef HMODULE lib_t;
#else
    #include <dlfcn.h>
    typedef void* lib_t;
#endif

namespace cdstar {

class Argument;
class NamedAssignment;
class StructType;

class Library {
public:
    Library(const std::string &name, bool includeHeaders = true);
    virtual ~Library();
    void AddStruct(const StructType &structType);
    void AddFunction(const std::vector<std::shared_ptr<Argument>> &input, 
                     const std::vector<std::shared_ptr<NamedAssignment>> &output,
                     const std::string &name);
    void CompileAndLoad();
    void* LoadFunction(const std::string &name);
private:
    std::fstream m_Stream;
    std::string m_Name;
    std::string m_CSourceFileName;
    std::string m_LibFileName;
    std::vector<std::string> m_FunctionNames;
    std::unordered_set<std::string> m_RegisteredStruct;
    lib_t m_LibHandle;
    bool m_Loaded;
};

}

#endif //CDSTAR_LIBRARY_H__