#include "library.h"
#include "expression.h"

namespace cdstar {

// http://stackoverflow.com/questions/2538103/how-to-call-a-function-from-a-shared-library
lib_t LoadLib(const char* lib) {
#ifdef _WIN32
    return ::LoadLibraryA(lib);
#else //_WIN32
    return ::dlopen(lib, RTLD_LAZY);
#endif //_WIN32
}

void UnloadLib(lib_t handle) {
#ifdef _WIN32
    ::FreeLibrary(handle);
#else //_WIN32
    ::dlclose(handle);
#endif //_WIN32
}

void* LoadFunction(lib_t handle, const char* name) {
# ifdef _WIN32
    return (void*)::GetProcAddress(handle, name);
# else //_WIN32
    return ::dlsym(handle, name);
# endif //_WIN32
}

Library::Library(const std::string &name, bool includeHeaders) : m_Name(name) {
    m_CSourceFileName = m_Name + ".c";    
    m_Stream.open(m_CSourceFileName.c_str(), std::fstream::out);  
    if (includeHeaders) {
        m_Stream << "#include <math.h>" << std::endl;
    }
    
#ifdef _WIN32    
    m_LibFileName = "lib" + m_Name + ".dll";
#else
    m_LibFileName = "lib" + m_Name + ".so";
#endif    
    m_Loaded = false;
}

Library::~Library() {    
    m_Stream.close();
    if (m_Loaded) {
        UnloadLib(m_LibHandle);
    }
}

void Library::AddStruct(const StructType &structType) {
    structType.EmitForwardDeclaration(m_Stream);
}

void Library::AddFunction(const std::vector<std::shared_ptr<Argument>> &input, 
                          const std::vector<std::shared_ptr<NamedAssignment>> &output,
                          const std::string &name) {
    EmitFunction(input, output, name, m_Stream);    
}

void Library::CompileAndLoad() {
    m_Stream.close();
#ifdef _WIN32    
    std::string cmd = "gcc -O3 -shared -o " + m_LibFileName + " " + m_CSourceFileName;
#else
    std::string cmd = "gcc -O3 -shared -fPIC -o " + m_LibFileName + " " + m_CSourceFileName;
#endif        
    std::system(cmd.c_str());
    m_LibHandle = LoadLib(m_LibFileName.c_str());    
    m_Loaded = true;
}

void* Library::LoadFunction(const std::string &name) {
    if (!m_Loaded) {
        throw std::runtime_error("Library Not Loaded");
    }
    return cdstar::LoadFunction(m_LibHandle, name.c_str());
}

} // namespace cdstar