#include "cfg.h"
#include "expression.h"
#include <unordered_set>

namespace cdstar {

void BuildCFGBlockRecursive(const Expression *expr,
                            std::shared_ptr<CFGBlock> &block,
                            std::unordered_set<const Expression*> &dfsSet) {
    if (!expr->UseTmpVar() && expr->Type() != ET_NAMED_ASSIGNMENT) {
        return;
    }
    if (dfsSet.find(expr) != dfsSet.end()) {
        return;
    }
    dfsSet.insert(expr);
    for (auto child : expr->Children()) {
        BuildCFGBlockRecursive(child.get(), block, dfsSet);
    }
    if (expr->Type() != ET_CONDEXPR) {
        block->exprs.push_back(expr);
    } else {
        std::shared_ptr<CFGBlock> trueBlock = std::make_shared<CFGBlock>();
        std::shared_ptr<CFGBlock> falseBlock = std::make_shared<CFGBlock>();
        std::shared_ptr<CFGBlock> newBlock = std::make_shared<CFGBlock>();
        block->nextBlocks = {trueBlock, falseBlock};
        block->cond = dynamic_cast<const Boolean*>(expr->Children()[0].get());
        block->condExprs.push_back(dynamic_cast<const CondExpr*>(expr));
        trueBlock->exprs = {expr};
        trueBlock->condId = 0;
        trueBlock->nextBlocks = {newBlock};
        trueBlock->prevBlocks = {block};
        falseBlock->exprs = {expr};
        falseBlock->condId = 1;
        falseBlock->nextBlocks = {newBlock};
        falseBlock->prevBlocks = {block};
        newBlock->prevBlocks = {trueBlock, falseBlock};
        block = newBlock;
    }
}

std::shared_ptr<CFGBlock> BuildCFGBlock(const std::vector<std::shared_ptr<NamedAssignment>> &exprs) {
    std::shared_ptr<CFGBlock> first = std::make_shared<CFGBlock>();
    first->isFirst = true;
    std::shared_ptr<CFGBlock> tmp = first;
    std::unordered_set<const Expression*> dfsSet;
    for (auto expr : exprs) {
        BuildCFGBlockRecursive(expr.get(), tmp, dfsSet);
    }
    return first;
}

bool Mergeable(std::shared_ptr<CFGBlock> ref,
               std::shared_ptr<CFGBlock> block,
               const std::vector<int> &nextId,
               const std::vector<std::shared_ptr<CFGBlock>> &blockList,
               const std::unordered_map<CFGBlock*, int> &blockMap,
               const std::set<std::pair<CFGBlock*, CFGBlock*>> &blockDependency) {
    if (block == ref) {
        return false;
    }
    if (ref->cond != block->cond) {
        return false;
    }
    return blockDependency.find({block.get(), ref.get()}) == blockDependency.end();
}

std::shared_ptr<CFGBlock> FindFirstMergeable(std::shared_ptr<CFGBlock> ref,
                                             const std::vector<int> &nextId,
                                             const std::vector<std::shared_ptr<CFGBlock>> &blockList,
                                             const std::unordered_map<CFGBlock*, int> &blockMap,
                                             const std::set<std::pair<CFGBlock*, CFGBlock*>> &blockDependency) {
    int refId = blockMap.at(ref.get());
    for (int i = nextId[refId]; i != -1; i = nextId[i]) {
        std::shared_ptr<CFGBlock> block = blockList[i];
        if (Mergeable(ref, block, nextId, blockList, blockMap, blockDependency)) {
            return block;
        }
    }
    return nullptr;
}

void BuildNextIdList(const std::shared_ptr<CFGBlock> block,
                     const std::unordered_map<CFGBlock*, int> &blockMap,
                     std::vector<int> &nextId,
                     std::stack<int> &stackId) {
    if (block->nextBlocks.size() == 1 && block->condId == 1) {
        int prevId = stackId.top();
        stackId.pop();
        nextId[prevId] = blockMap.at(block->nextBlocks[0].get());
        BuildNextIdList(block->nextBlocks[0], blockMap, nextId, stackId);
    }
    if (block->nextBlocks.size() == 2) {
        stackId.push(blockMap.at(block.get()));
        BuildNextIdList(block->nextBlocks[0], blockMap, nextId, stackId);
        BuildNextIdList(block->nextBlocks[1], blockMap, nextId, stackId);
    }
}

// block A depends on B iff: A depends on conditions stored in B
void BuildBlockDependency(const std::vector<std::shared_ptr<CFGBlock>> &blockList,
                          const std::unordered_map<CFGBlock*, int> &blockMap,
                          const std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> &exprMap,
                          std::set<std::pair<CFGBlock*, CFGBlock*>> &blockDependency) {
    for (auto block : blockList) {
        for (auto condExpr : block->condExprs) {
            std::unordered_set<const Expression*> dfsSet;
            std::stack<const Expression*> dfsStack;
            dfsStack.push(condExpr);
            while (!dfsStack.empty()) {
                const Expression *expr = dfsStack.top();
                dfsStack.pop();
                if (dfsSet.find(expr) != dfsSet.end()) {
                    continue;
                }
                dfsSet.insert(expr);

                if (expr->Type() == ET_CONDEXPR && expr != condExpr) {
                    blockDependency.insert({block.get(), exprMap.at(expr).get()});
                }

                for (auto child : expr->Children()) {
                    dfsStack.push(child.get());
                }
            }
        }
    }
}

void BuildExprBlockMap(const std::vector<std::shared_ptr<CFGBlock>> &blockList,
                       std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> &exprMap) {
    for (auto block : blockList) {
        for (auto expr : block->exprs) {
            if (expr->Type() != ET_CONDEXPR) {
                exprMap[expr] = block;
            }
        }
        for (auto condExpr : block->condExprs) {
            exprMap[condExpr] = block;
        }
    }
}

void BuildCFGBlockList(const std::shared_ptr<CFGBlock> block,
                       std::vector<std::shared_ptr<CFGBlock>> &blockList,
                       std::unordered_map<CFGBlock*, int> &blockMap) {
    blockList.push_back(block);
    blockMap[block.get()] = blockList.size() - 1;
    if (block->nextBlocks.size() == 1 && block->condId == 1) {
        BuildCFGBlockList(block->nextBlocks[0], blockList, blockMap);
    }
    if (block->nextBlocks.size() == 2) {
        BuildCFGBlockList(block->nextBlocks[0], blockList, blockMap);
        BuildCFGBlockList(block->nextBlocks[1], blockList, blockMap);
    }
}

void PushExpr(const Expression *expr,
              const std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> &exprMap,
              const std::unordered_map<CFGBlock*, int> &blockMap,
              std::shared_ptr<CFGBlock> toPush,
              std::unordered_set<const Expression*> &dfsSet,
              std::unordered_set<const Expression*> &insertedSet,
              int toMergeBlockId,
              int mergeBlockId) {
    if (dfsSet.find(expr) != dfsSet.end() ||
        insertedSet.find(expr) != insertedSet.end()) {
        return;
    }

    dfsSet.insert(expr);
    
    for (auto child : expr->Children()) {
        PushExpr(child.get(), exprMap, blockMap, toPush, dfsSet, insertedSet, toMergeBlockId, mergeBlockId);
    }

    if (expr->Type() != ET_CONDEXPR && expr->UseTmpVar()) {
        std::shared_ptr<CFGBlock> exprBlock = exprMap.at(expr);
        int exprBlockId = blockMap.at(exprBlock.get());
        if (exprBlockId > toMergeBlockId && exprBlockId <= mergeBlockId) {
            insertedSet.insert(expr);
            toPush->exprs.push_back(expr);
            exprBlock->exprs.remove(expr);
        }
    }
}

void Merge(std::shared_ptr<CFGBlock> toMerge,
           std::shared_ptr<CFGBlock> merge,
           const std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> &exprMap,
           const std::vector<std::shared_ptr<CFGBlock>> &blockList,
           const std::unordered_map<CFGBlock*, int> &blockMap,
           const std::vector<int> &nextId) {
    int toMergeBlockId = blockMap.at(toMerge.get());
    int mergeBlockId = blockMap.at(merge.get());
    std::unordered_set<const Expression*> insertedSet;
    // move all the expressions depended by the expressions in the true/false clauses, 
    // in the the toMerge block
    for (int i = mergeBlockId + 1; i < nextId[mergeBlockId]; i++) {
        std::shared_ptr<CFGBlock> curBlock = blockList[i];
        for (auto expr : curBlock->exprs) {
            std::unordered_set<const Expression*> dfsSet;
            if (expr->Type() != ET_CONDEXPR) {
                for (auto child : expr->Children()) {
                    PushExpr(child.get(), exprMap, blockMap, toMerge, dfsSet, insertedSet, toMergeBlockId, mergeBlockId);
                }
            } else {
                PushExpr(expr->Children()[curBlock->condId+1].get(), 
                    exprMap, blockMap, toMerge, dfsSet, insertedSet, toMergeBlockId, mergeBlockId);
            }
        }
    }

    std::shared_ptr<CFGBlock> toMergeLastBlock = blockList[nextId[toMergeBlockId]];
    std::shared_ptr<CFGBlock> mergeLastBlock = blockList[nextId[mergeBlockId]];
    std::shared_ptr<CFGBlock> toMergeTrueBlock = toMergeLastBlock->prevBlocks[0].lock();
    for (auto expr : merge->nextBlocks[0]->exprs) {
        toMergeTrueBlock->exprs.push_back(expr);
    }
    if (merge->nextBlocks[0]->nextBlocks.size() > 1) {
        toMergeTrueBlock->nextBlocks = merge->nextBlocks[0]->nextBlocks;
        toMergeTrueBlock->condExprs = merge->nextBlocks[0]->condExprs;
        std::shared_ptr<CFGBlock> mergeTrueBlock = mergeLastBlock->prevBlocks[0].lock();
        mergeTrueBlock->nextBlocks = {std::shared_ptr<CFGBlock>(toMergeLastBlock)};
        toMergeLastBlock->prevBlocks[0] = mergeLastBlock->prevBlocks[0];
    }
    std::shared_ptr<CFGBlock> toMergeFalseBlock = toMergeLastBlock->prevBlocks[1].lock();
    for (auto expr : merge->nextBlocks[1]->exprs) {
        toMergeFalseBlock->exprs.push_back(expr);
    }
    if (merge->nextBlocks[1]->nextBlocks.size() > 1) {
        toMergeFalseBlock->nextBlocks = merge->nextBlocks[1]->nextBlocks;
        toMergeFalseBlock->condExprs = merge->nextBlocks[1]->condExprs;
        std::shared_ptr<CFGBlock> mergeFalseBlock = mergeLastBlock->prevBlocks[1].lock();
        mergeFalseBlock->nextBlocks = {std::shared_ptr<CFGBlock>(toMergeLastBlock)};
        toMergeLastBlock->prevBlocks[1] = mergeLastBlock->prevBlocks[1];
    }
    toMerge->condExprs.insert(toMerge->condExprs.end(), merge->condExprs.begin(), merge->condExprs.end());
    for (auto expr : mergeLastBlock->exprs) {
        merge->exprs.push_back(expr);
    }
    merge->nextBlocks = mergeLastBlock->nextBlocks;
    if (mergeLastBlock->nextBlocks.size() > 0) {
        mergeLastBlock->prevBlocks = {std::shared_ptr<CFGBlock>(merge)};
    }
    merge->cond = mergeLastBlock->cond;
    merge->condExprs = mergeLastBlock->condExprs;
}

struct CFGBlockTreeNode {
    std::weak_ptr<CFGBlockTreeNode> parent;
    std::vector<std::shared_ptr<CFGBlockTreeNode>> children;
    std::vector<std::shared_ptr<CFGBlock>> blocks;
};

std::shared_ptr<CFGBlockTreeNode> BuildCFGBlockTree(const std::vector<std::shared_ptr<CFGBlock>> &blockList,
                                                    std::unordered_map<CFGBlock*, std::shared_ptr<CFGBlockTreeNode>> &blockNodeMap) {
    std::shared_ptr<CFGBlockTreeNode> root = std::make_shared<CFGBlockTreeNode>();
    std::stack<std::shared_ptr<CFGBlockTreeNode>> stackNodes;
    stackNodes.push(root);
    for (auto block : blockList) {
        auto curNode = stackNodes.top();
        curNode->blocks.push_back(block);
        blockNodeMap[block.get()] = curNode;
        if (block->nextBlocks.size() == 2) {
            std::shared_ptr<CFGBlockTreeNode> newNode = std::make_shared<CFGBlockTreeNode>();
            curNode->children.push_back(newNode);
            newNode->parent = curNode;
            stackNodes.push(newNode);
        }
        if (block->nextBlocks.size() == 1 && block->condId == 0) {
            stackNodes.pop();
            std::shared_ptr<CFGBlockTreeNode> newNode = std::make_shared<CFGBlockTreeNode>();
            stackNodes.top()->children.push_back(newNode);
            newNode->parent = stackNodes.top();
            stackNodes.push(newNode);
        }
        if (block->nextBlocks.size() == 1 && block->condId == 1) {
            stackNodes.pop();
        }
    }
    return root;
}

std::shared_ptr<CFGBlockTreeNode> LowestCommonAncestor(std::shared_ptr<CFGBlockTreeNode> a,
                                                       std::shared_ptr<CFGBlockTreeNode> b) {
    if (a == b) {
        return a;
    }
    std::vector<std::shared_ptr<CFGBlockTreeNode>> aParents = {a};
    std::vector<std::shared_ptr<CFGBlockTreeNode>> bParents = {b};
    for (auto p = a->parent.lock(); p != nullptr; p = p->parent.lock()) {
        aParents.push_back(p);
    }
    for (auto p = b->parent.lock(); p != nullptr; p = p->parent.lock()) {
        bParents.push_back(p);
    }
    for (auto aP : aParents) {
        for (auto bP : bParents) {
            if (aP == bP) {
                return aP;
            }
        }
    }
    return nullptr;
}

void BuildCondBlockDep(const std::vector<std::shared_ptr<CFGBlock>> &blockList,
                       std::unordered_map<CFGBlock*, std::shared_ptr<CFGBlockTreeNode>> &blockNodeMap,
                       std::unordered_map<const Expression*, std::shared_ptr<CFGBlockTreeNode>> &condNodeDep) {
    std::shared_ptr<CFGBlockTreeNode> root = BuildCFGBlockTree(blockList, blockNodeMap);
    std::unordered_set<const Expression*> exprSet;
    std::unordered_multimap<const Expression*, std::shared_ptr<CFGBlock>> condBlockDep;
    for (auto block : blockList) {
        std::unordered_set<const Expression*> visited;
        for (auto expr : block->exprs) {
            if (expr->Type() != ET_CONDEXPR) {
                for (auto child : expr->Children()) {
                    if (visited.find(child.get()) != visited.end()) {
                        continue;
                    }
                    condBlockDep.insert({child.get(), block});
                    exprSet.insert(child.get());
                    visited.insert(child.get());
                }
            } else {
                auto child = expr->Children()[block->condId+1];
                if (visited.find(child.get()) != visited.end()) {
                    continue;
                }
                condBlockDep.insert({child.get(), block});
                exprSet.insert(child.get());
                visited.insert(child.get());
            }
        }
    }

    for (auto expr : exprSet) {
        std::unordered_set<std::shared_ptr<CFGBlockTreeNode>> nodes;
        for (auto it = condBlockDep.equal_range(expr); it.first != it.second; it.first++) {
            auto keyValue = *(it.first);
            nodes.insert(blockNodeMap[keyValue.second.get()]);
            for (int i = 0; i < blockList.size(); i++) {
                if (blockList[i] == keyValue.second) {
                    break;
                }
            }
        }

        if (nodes.size() == 1) {
            condNodeDep[expr] = *(nodes.begin());
        } else if (nodes.size() >= 2) {
            std::shared_ptr<CFGBlockTreeNode> n = *(nodes.begin());
            for (auto node : nodes) {
                n = LowestCommonAncestor(n, node);
            }
            condNodeDep[expr] = n;
        }
    }
}

bool PullExpr(const std::vector<std::shared_ptr<CFGBlock>> &blockList,
              const std::unordered_map<CFGBlock*, int> &blockMap,
              const std::vector<int> &nextId,
              const std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> &exprMap) {
    std::unordered_map<CFGBlock*, std::shared_ptr<CFGBlockTreeNode>> blockNodeMap;
    std::unordered_map<const Expression*, std::shared_ptr<CFGBlockTreeNode>> condNodeDep;
    BuildCondBlockDep(blockList, blockNodeMap, condNodeDep);
    bool changed = false;
    std::set<std::pair<const Expression*, const CFGBlockTreeNode*>> pulledExpr;
    for (int blockId = 0; blockId < (int)blockList.size(); blockId++) {
        auto block = blockList[blockId];
        if (block->condId == -1) {
            continue;
        }
        for (auto expr : block->exprs) {
            for (auto child : expr->Children()) {
                if (expr->Type() == ET_CONDEXPR && child != expr->Children()[block->condId+1]) {
                    continue;
                }
                auto exprMapIt = exprMap.find(child.get());
                if (exprMapIt == exprMap.end()) {
                    continue;
                }
                auto exprBlock = exprMapIt->second;
                if (exprBlock == block) {
                    continue;
                }
                auto condNodeIt = condNodeDep.find(child.get());
                if (condNodeIt == condNodeDep.end()) {
                    continue;
                }
                auto condNode = condNodeIt->second;
                if (child->Type() == ET_CONDEXPR) {
                    auto nodeAncestor = condNode;
                    for (auto condExpr : exprBlock->condExprs) {
                        auto condNodeIt = condNodeDep.find(condExpr);
                        if (condNodeIt == condNodeDep.end()) {
                            continue;
                        }
                        nodeAncestor = LowestCommonAncestor(nodeAncestor, condNodeIt->second);
                    }
                    if (nodeAncestor == nullptr) {
                        continue;
                    }
                    if (blockNodeMap.find(exprBlock.get())->second == nodeAncestor) {
                        continue;
                    }
                    auto moveInBlock = nodeAncestor->blocks[0];

                    int exprBlockNextId = nextId[blockMap.find(exprBlock.get())->second];
                    auto exprBlockNext = blockList[exprBlockNextId];
                    std::shared_ptr<CFGBlock> newBlock = std::make_shared<CFGBlock>();
                    newBlock->nextBlocks = exprBlock->nextBlocks;
                    newBlock->prevBlocks = moveInBlock->prevBlocks;
                    newBlock->condId = moveInBlock->condId;
                    newBlock->cond = exprBlock->cond;
                    newBlock->condExprs = exprBlock->condExprs;
                    newBlock->nextBlocks[0]->prevBlocks = {newBlock};
                    newBlock->nextBlocks[1]->prevBlocks = {newBlock};
                    
                    if (moveInBlock->prevBlocks.size() == 1) {
                        moveInBlock->prevBlocks[0].lock()->nextBlocks[moveInBlock->condId] = newBlock;
                    } else if (moveInBlock->prevBlocks.size() == 2) {
                        moveInBlock->prevBlocks[0].lock()->nextBlocks = {newBlock};
                        moveInBlock->prevBlocks[1].lock()->nextBlocks = {newBlock};
                    }
                    moveInBlock->prevBlocks = exprBlockNext->prevBlocks;

                    exprBlock->exprs.insert(exprBlock->exprs.end(), exprBlockNext->exprs.begin(), exprBlockNext->exprs.end());
                    exprBlock->nextBlocks = exprBlockNext->nextBlocks;
                    exprBlock->cond = exprBlockNext->cond;
                    exprBlock->condExprs = exprBlockNext->condExprs;
                    for (auto nextBlock : exprBlockNext->nextBlocks) {
                        nextBlock->prevBlocks = {exprBlock};
                    }
                    
                    exprBlockNext->exprs.clear();
                    exprBlockNext->prevBlocks[0].lock()->nextBlocks[0] = moveInBlock;
                    exprBlockNext->prevBlocks[1].lock()->nextBlocks[0] = moveInBlock;
                    return true;
                } else {
                    if (blockNodeMap.find(exprBlock.get())->second == condNode) {
                        continue;
                    }
                    if (pulledExpr.find({child.get(), condNode.get()}) != pulledExpr.end()) {
                        continue;
                    }
                    auto moveInBlock = condNode->blocks[0];
                    pulledExpr.insert({child.get(), condNode.get()});
                    exprBlock->exprs.remove(child.get());
                    moveInBlock->exprs.push_front(child.get());
                    changed = true;
                }
            }
        }
    }
    return changed;
}

#if 0
bool OptimizeNestedCond(const std::vector<std::shared_ptr<CFGBlock>> &blockList) {
    std::deque<std::pair<const Boolean*,int>> condStack;
    for (int blockId = 0; blockId < (int)blockList.size(); blockId++) {
        auto block = blockList[blockId];
        if (block->nextBlocks.size() == 2) {
            for (auto cond : condStack) {
                if (cond.first == block->cond) {
                    block->exprs.insert(block->exprs.end(), 
                        block->nextBlocks[cond.second]->exprs.begin(),
                        block->nextBlocks[cond.second]->exprs.end());
                    if (block->nextBlocks[cond.second]->nextBlocks.size() == 1) {
                        block->exprs.insert(block->exprs.end(),
                            block->nextBlocks[cond.second]->nextBlocks[0]->exprs.begin(),
                            block->nextBlocks[cond.second]->nextBlocks[0]->exprs.end());
                        block->nextBlocks = block->nextBlocks[cond.second]->nextBlocks[0]->nextBlocks;
                        block->cond = block->nextBlocks[cond.second]->nextBlocks[0]->cond;
                        block->condExprs = block->nextBlocks[cond.second]->nextBlocks[0]->condExprs;
                    } else if (block->nextBlocks[cond.second]->nextBlocks.size() == 2) {
                        block->nextBlocks = block->nextBlocks[cond.second]->nextBlocks;
                        block->cond = block->nextBlocks[cond.second]->cond;
                        block->condExprs = block->nextBlocks[cond.second]->condExprs;
                    }
                    for (auto nextBlock : block->nextBlocks) {
                        nextBlock->prevBlocks = {block};
                    }
                    return true;
                }
            }
            condStack.push_back({block->cond, 0});
        }
        if (block->nextBlocks.size() == 1 && block->condId == 0) {
            const Boolean* cond = condStack.back().first;
            condStack.pop_back();
            condStack.push_back({cond, 1});
        }
        if (block->nextBlocks.size() == 1 && block->condId == 1) {
            condStack.pop_back();
        }
    }
    return false;
}
#endif

void Optimize(std::shared_ptr<CFGBlock> block) {
    bool changed = false;
    do {
        changed = false;
        std::vector<std::shared_ptr<CFGBlock>> blockList;
        std::unordered_map<CFGBlock*, int> blockMap;
        std::unordered_map<const Expression*, std::shared_ptr<CFGBlock>> exprMap;
        std::set<std::pair<CFGBlock*, CFGBlock*>> blockDependency;
        BuildCFGBlockList(block, blockList, blockMap);
        BuildExprBlockMap(blockList, exprMap);
        BuildBlockDependency(blockList, blockMap, exprMap, blockDependency);
        std::vector<int> nextId;
        nextId.resize(blockList.size(), -1);
        std::stack<int> stackId;
        BuildNextIdList(block, blockMap, nextId, stackId);
        if (!changed) {
            changed = PullExpr(blockList, blockMap, nextId, exprMap);
        }
        if (!changed) {
            for (auto toMerge : blockList) {
                std::shared_ptr<CFGBlock> merge = FindFirstMergeable(toMerge, nextId, blockList, blockMap, blockDependency);
                if (merge != nullptr) {
                    Merge(toMerge, merge, exprMap, blockList, blockMap, nextId);
                    changed = true;
                    break;
                }
            }
        }
    } while(changed);
}

void BuildExprMap(const std::shared_ptr<CFGBlock> block,
                  std::unordered_map<const Expression*, ExprInfo> &exprMap,
                  int &tmpId,
                  int &boolId,
                  int depth = 1) {
    for (auto expr : block->exprs) {
        if (exprMap.find(expr) != exprMap.end() || expr->Type() == ET_NAMED_ASSIGNMENT) {
            continue;
        }
        ExprInfo exprInfo;
        if (expr->Type() != ET_BOOLEAN) {
            exprInfo.index = tmpId++;
        } else {
            exprInfo.index = boolId++;
        }
        exprMap[expr] = exprInfo;
    }
    if (block->nextBlocks.size() == 1 && block->condId == 1) {
        BuildExprMap(block->nextBlocks[0], exprMap, tmpId, boolId, depth - 1);
    }
    if (block->nextBlocks.size() == 2) {
        BuildExprMap(block->nextBlocks[0], exprMap, tmpId, boolId, depth + 1);
        BuildExprMap(block->nextBlocks[1], exprMap, tmpId, boolId, depth + 1);
    }
}

void EmitCFGBlock(const std::shared_ptr<CFGBlock> block,
                  std::unordered_map<const Expression*, ExprInfo> &exprMap,
                  std::ostream &os,
                  std::deque<std::pair<std::shared_ptr<CFGBlock>, bool>> &blockStack,
                  int depth = 1) {
    bool nested = false;
    bool hideEmit = false;
    if (block->condId >= 0 && blockStack.size() >= 2) {
        auto topBlock = blockStack.back().first;
        for (int i = 0; i < (int)blockStack.size() - 1; i++) {
            auto b = blockStack[i];
            if (b.first->cond == topBlock->cond) {
                nested = true;
                hideEmit = (b.second && block->condId == 1) || (!b.second && block->condId == 0);
                break;
            }
        }
    }
    if (!hideEmit) {
        for (auto expr : block->exprs) {
            PrintTab(nested ? depth - 1 : depth, os);
            expr->Emit(block.get(), exprMap, os);
        }
    }
    if (block->nextBlocks.size() == 1 && block->condId == 0) {
        if (!nested) {
            PrintTab(depth - 1, os) << "} else { //" << blockStack.back().first->cond->GetEmitName(exprMap) << std::endl;
        }
        auto topBlock = blockStack.back().first;
        blockStack.pop_back();
        blockStack.push_back({topBlock, false});
    }
    if (block->nextBlocks.size() == 1 && block->condId == 1) {
        if (!nested) {
            PrintTab(depth - 1, os) << "} //" << blockStack.back().first->cond->GetEmitName(exprMap) << std::endl;
        }
        blockStack.pop_back();
        EmitCFGBlock(block->nextBlocks[0], exprMap, os, blockStack, depth - 1);
    }
    if (block->nextBlocks.size() == 2) {
        bool nested = false;
        for (auto b : blockStack) {
            if (b.first->cond == block->cond) {
                nested = true;
                break;
            }
        }
        if (!nested) {
            std::string condName = block->cond->GetEmitName(exprMap);
            PrintTab(depth, os) << "if (" << condName << ") {" << std::endl;
        }
        blockStack.push_back({block, true});
        EmitCFGBlock(block->nextBlocks[0], exprMap, os, blockStack, depth + 1);
        EmitCFGBlock(block->nextBlocks[1], exprMap, os, blockStack, depth + 1);
    }
}

void Emit(const std::shared_ptr<CFGBlock> block, std::ostream &os) {
    std::unordered_map<const Expression*, ExprInfo> exprMap;
    int tmpId = 0;
    int boolId = 0;
    BuildExprMap(block, exprMap, tmpId, boolId);
    os << "\tdouble _t[" << tmpId << "];" << std::endl;
    std::deque<std::pair<std::shared_ptr<CFGBlock>,bool>> blockStack;
    EmitCFGBlock(block, exprMap, os, blockStack);
}

}