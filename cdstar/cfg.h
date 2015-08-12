#ifndef CDSTAR_CFG_H__
#define CDSTAR_CFG_H__

#include <list>
#include <vector>
#include <ostream>

namespace cdstar {

class Expression;
class Boolean;
class NamedAssignment;
class CondExpr;
    
struct CFGBlock {
    CFGBlock() {
        isFirst = false;
        condId = -1;
        cond = nullptr;
    }

    bool isFirst;
    std::list<const Expression*> exprs;
    std::vector<std::shared_ptr<CFGBlock>> nextBlocks;
    std::vector<std::weak_ptr<CFGBlock>> prevBlocks;
    int condId;
    const Boolean *cond;
    std::vector<const CondExpr*> condExprs;
};

struct ExprInfo {
    int index;
};

std::shared_ptr<CFGBlock> BuildCFGBlock(const std::vector<std::shared_ptr<NamedAssignment>> &exprs);
void Optimize(std::shared_ptr<CFGBlock> block);
void Emit(const std::shared_ptr<CFGBlock> block, std::ostream &os);

}

#endif // CDSTAR_CFG_H__