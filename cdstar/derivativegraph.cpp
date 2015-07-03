#include "derivativegraph.h"
#include "expression.h"

#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <functional>
#include <set>

namespace cdstar {

void PrintGraph(const DervGraph& graph) {
    auto indexMap = boost::get(boost::vertex_index, graph);
    for(auto it = boost::vertices(graph); it.first != it.second; it.first++) {
        auto v = *it.first;
        std::cerr << indexMap[v] << ":" << std::endl;
        for(auto it = boost::out_edges(v, graph); it.first != it.second; it.first++) {
            auto e = *it.first;
            std::cerr << indexMap[boost::target(e, graph)] << std::endl;
            std::shared_ptr<Expression> expr = boost::get(edge_property_expression_t(), graph, e);
            expr->Print();
        }
    }
}

class DominatorsFinder {
public:
    DominatorsFinder(const DervGraph&       graph, 
                     const DervGraphVertex& start, 
                     const bool             post);
    std::vector<int> GetDominatorsMap() const {
        return m_Dominators;
    }

private:
    void BuildPostOrderMap(const DervGraph&       graph, 
                           const DervGraphVertex& vert, 
                           const bool             post,
                           int&                   number);
    int Intersect(int b1, int b2);

    std::vector<int> m_PostOrderMap;
    std::vector<int> m_InvPostOrderMap;
    std::vector<int> m_Dominators;
};

void DominatorsFinder::BuildPostOrderMap(const DervGraph&       graph, 
                                         const DervGraphVertex& vert, 
                                         const bool             post,
                                         int&                   number) {    
    if(m_PostOrderMap[vert] != -1) {
        return;
    }
    if(post) {
        for(auto it = boost::in_edges(vert, graph); it.first != it.second; it.first++) {
            BuildPostOrderMap(graph, boost::source(*it.first, graph), post, number);
        }
    } else {
        for(auto it = boost::out_edges(vert, graph); it.first != it.second; it.first++) {
            BuildPostOrderMap(graph, boost::target(*it.first, graph), post, number);
        }
    }
    m_InvPostOrderMap[number] = vert;
    m_PostOrderMap[vert] = number;
    number++;
}

DominatorsFinder::DominatorsFinder(const DervGraph&       graph, 
                                   const DervGraphVertex& start, 
                                   const bool             post) {
    // compute post order numbering
    m_PostOrderMap.resize(boost::num_vertices(graph), -1);
    m_InvPostOrderMap.resize(boost::num_vertices(graph), -1);
    int number = 0;
    BuildPostOrderMap(graph, start, post, number);    
    // initialize the dominators array
    m_Dominators.resize(boost::num_vertices(graph), -1);
    m_Dominators[start] = start;
    int startPostIndex = m_PostOrderMap[start];
    bool changed = true;
    while(changed) {
        changed = false;
        // for all nodes, b, in reverse postorder (except start_node)
        for(int postId = (int)m_Dominators.size() - 1; postId >= 0; postId--) {
            if(postId == startPostIndex) {
                continue;
            }
            if(m_InvPostOrderMap[postId] == -1) { // not connected to starting point
                continue;
            }
            auto b = boost::vertex(m_InvPostOrderMap[postId], graph);
            int newIDom = -1;
            if(post) {
                for(auto it = boost::out_edges(b, graph); it.first != it.second; it.first++) {
                    auto e = *it.first;
                    auto p = boost::target(e, graph);
                    if(m_Dominators[p] != -1) { // processed
                        if(newIDom == -1) {
                            newIDom = p;
                        } else {
                            newIDom = Intersect(p, newIDom);
                        }
                    }
                }
            } else {
                for(auto it = boost::in_edges(b, graph); it.first != it.second; it.first++) {
                    auto e = *it.first;
                    auto p = boost::source(e, graph);
                    if(m_Dominators[p] != -1) { // processed
                        if(newIDom == -1) {
                            newIDom = p;
                        } else {
                            newIDom = Intersect(p, newIDom);
                        }
                    }
                }
            }

            if(m_Dominators[b] != newIDom) {
                m_Dominators[b] = newIDom;
                changed = true;
            }
        }
    }
}

int DominatorsFinder::Intersect(int b1, int b2) {
    int finger1 = b1;
    int finger2 = b2;
    while(finger1 != finger2) {
        while(m_PostOrderMap[finger1] < m_PostOrderMap[finger2]) {
            finger1 = m_Dominators[finger1];
        }
        while(m_PostOrderMap[finger2] < m_PostOrderMap[finger1]) {
            finger2 = m_Dominators[finger2];
        }
    }
    return finger1;
}

std::pair<bool, std::shared_ptr<Expression>> 
        FactorExpr(const DervGraph                &graph,
                   const DervGraphVertex          &vert0,
                   const DervGraphVertex          &vert1,
                   std::set<std::pair<int, int> > &edges) {
    if(vert0 == vert1) {
        return std::pair<bool, std::shared_ptr<Expression> >(true, std::make_shared<Constant>(1.0));
    }    
    bool primary = false;
    std::shared_ptr<Expression> expr = std::make_shared<Constant>(0.0);
    for(auto it = boost::out_edges(vert0, graph); it.first != it.second; it.first++) {
        auto e = *it.first;
        std::shared_ptr<Expression> edgeExpr = boost::get(edge_property_expression_t(), graph, e);
        bool primaryEdge;
        std::shared_ptr<Expression> factorExpr;
        std::tie(primaryEdge, factorExpr) = FactorExpr(graph, boost::target(e, graph), vert1, edges);
        if(primaryEdge) {
            edges.insert(std::make_pair<int, int>(vert0, boost::target(e, graph)));
            primary = true;
            expr = expr + edgeExpr * factorExpr;
        }
    }
    return std::pair<bool, std::shared_ptr<Expression> >(primary, expr);
}

DerivativeGraph::DerivativeGraph(const std::vector<ExprPtrPair>& dervExprs) {
    std::vector<std::pair<DervGraphVertex, DervGraphVertex>> srcTgtPairs(dervExprs.size());
    for(size_t i = 0; i < dervExprs.size(); i++) {
        srcTgtPairs[i].first = BuildGraph(dervExprs[i].first);
    }    
    for(size_t i = 0; i < dervExprs.size(); i++) {
        bool found = false;
        for(auto it = boost::vertices(m_Graph); it.first != it.second; it.first++) {
            auto expr = boost::get(vertex_property_expression_t(), m_Graph, *it.first);
            if(expr == dervExprs[i].second) {
                srcTgtPairs[i].second = *it.first;
                found = true;
                break;
            }
        }
        if(!found) {
            auto vertex = boost::add_vertex(m_Graph);
            boost::put(vertex_property_expression_t(), m_Graph, vertex, dervExprs[i].second);
            srcTgtPairs[i].second = vertex;
            auto edge = boost::add_edge(srcTgtPairs[i].first, srcTgtPairs[i].second, m_Graph).first;
            boost::put(edge_property_expression_t(), m_Graph, edge, std::make_shared<Constant>(0.0));
        }
    }

    std::vector<DerivativeSubgraph> subgraphs;
    for(auto srcTgtPair : srcTgtPairs) {
        const auto& source = srcTgtPair.first;
        const auto& target = srcTgtPair.second;
        DervGraph subgraph(boost::num_vertices(m_Graph));
        std::unordered_map<DervGraphVertex, bool> vertSet;
        BuildSubgraph(target, source, subgraph, vertSet);
        DerivativeSubgraph dSubgrpah = { subgraph, source, target };
        subgraphs.push_back(dSubgrpah);
    }            
    FactorSubgraphs(subgraphs);             
    FactorCommonSubproduct(subgraphs, m_Derivatives);          
}

bool DerivativeGraph::GetFactorSubgraph(std::vector<DerivativeSubgraph> &subgraphs,
                                        DervGraphVertex                 &vert0,
                                        DervGraphVertex                 &vert1,
                                        std::vector<int>                &subgraphIdList,
                                        std::vector<bool>               &isDomList) {
    struct FactorGraphStruct {
        FactorGraphStruct() {}
        FactorGraphStruct(const DervGraphVertex& vert0,
                          const DervGraphVertex& vert1,
                          const int              subgraphId,
                          const bool             isDom)
            : vert0(vert0), vert1(vert1) {
            subgraphIdList.push_back(subgraphId);            
            isDomList.push_back(isDom);
        }

        DervGraphVertex vert0;
        DervGraphVertex vert1;
        std::vector<int> subgraphIdList;
        std::vector<bool> isDomList;
    };
    int highestCount = -1;
    std::pair<int, int> highestPair;
    std::map<std::pair<int, int>, FactorGraphStruct> factorGraphMap;
    for(int subgraphId = 0; subgraphId < (int)subgraphs.size(); subgraphId++) {
        DerivativeSubgraph& dSubgraph = subgraphs[subgraphId];
        const DervGraph& subgraph = dSubgraph.subgraph;        
        const DervGraphVertex& target = dSubgraph.target;
        const DervGraphVertex& source = dSubgraph.source;
        std::vector<int>& doms = dSubgraph.doms;
        std::vector<int>& postDoms = dSubgraph.postDoms;
        DominatorsFinder domFinder(subgraph, source, false);
        doms = domFinder.GetDominatorsMap();
        DominatorsFinder postDomFinder(subgraph, target, true);
        postDoms = postDomFinder.GetDominatorsMap();                
        std::set<std::pair<int, int>> subgraphSet;
        for(int vertId = 0; vertId < (int)doms.size(); vertId++) {
            auto v0 = boost::vertex(vertId, subgraph);
            if(boost::in_degree(v0, subgraph) <= 1 || doms[vertId] == -1) {
                continue;
            }
            for(int domId = doms[vertId]; domId != doms[domId]; domId = doms[domId]) {
                auto v1 = boost::vertex(domId, subgraph);
                if(boost::out_degree(v1, subgraph) > 1) {                    
                    std::pair<int, int> vPair = {v1, v0};
                    if (subgraphSet.find(vPair) != subgraphSet.end()) {
                        continue;
                    }
                    subgraphSet.insert(vPair);
                    if(factorGraphMap.find(vPair) == factorGraphMap.end()) {                        
                        factorGraphMap[vPair] = FactorGraphStruct(v1, v0, subgraphId, true);
                    } else {
                        factorGraphMap[vPair].subgraphIdList.push_back(subgraphId);
                        factorGraphMap[vPair].isDomList.push_back(true);
                    }
                    if((int)factorGraphMap[vPair].subgraphIdList.size() > highestCount) {
                        highestCount = factorGraphMap[vPair].subgraphIdList.size();
                        highestPair = vPair;
                    }
                }
            }
        }

        for(int vertId = 0; vertId < (int)postDoms.size(); vertId++) {
            auto v0 = boost::vertex(vertId, subgraph);
            if(boost::out_degree(v0, subgraph) <= 1 || postDoms[vertId] == -1) {
                continue;
            }
            for(int domId = postDoms[vertId]; domId != postDoms[domId]; domId = postDoms[domId]) {
                auto v1 = boost::vertex(domId, subgraph);
                if(boost::in_degree(v1, subgraph) > 1) {
                    std::pair<int, int> vPair = {v0, v1};
                    if (subgraphSet.find(vPair) != subgraphSet.end()) {
                        continue;
                    }
                    subgraphSet.insert(vPair);
                    if(factorGraphMap.find(vPair) == factorGraphMap.end()) {
                        factorGraphMap[vPair] = FactorGraphStruct(v0, v1, subgraphId, false);
                    } else {
                        factorGraphMap[vPair].subgraphIdList.push_back(subgraphId);
                        factorGraphMap[vPair].isDomList.push_back(false);
                    }
                    if((int)factorGraphMap[vPair].subgraphIdList.size() > highestCount) {
                        highestCount = factorGraphMap[vPair].subgraphIdList.size();
                        highestPair = vPair;
                    }
                }
            }
        }
    }

    if(highestCount == -1) {
        return false;
    }

    FactorGraphStruct &fg = factorGraphMap[highestPair];
    vert0 = fg.vert0;
    vert1 = fg.vert1;
    subgraphIdList = fg.subgraphIdList;
    isDomList = fg.isDomList;

    return true;
}

DervGraphVertex DerivativeGraph::BuildGraph(const std::shared_ptr<Expression> node) {
    if(m_Vertices.find(node.get()) != m_Vertices.end()) {
        return m_Vertices[node.get()];
    }

    DervGraphVertex vertex = m_Vertices[node.get()] = boost::add_vertex(m_Graph);
    boost::put(vertex_property_expression_t(), m_Graph, vertex, node);
    
    std::unordered_map<DervGraphVertex, int> vertexRec;
    auto children = node->Children();
    auto dervs = node->Dervs();
    for(size_t i = 0; i < children.size(); i++) {
        auto child = children[i];
        auto derv = dervs[i];
        if (derv->Type() == ET_CONSTANT && derv->GetConstant() == 0.0) {
            continue;
        }
        DervGraphVertex childVertex = BuildGraph(child);
        vertexRec[childVertex] = vertexRec[childVertex] + 1;
        if (vertexRec[childVertex] == 1) {
            boost::add_edge(vertex, childVertex, derv, m_Graph);
        } else {
            auto edgeRet = boost::edge(vertex, childVertex, m_Graph);
            std::shared_ptr<Expression> expr = boost::get(edge_property_expression_t(), m_Graph, edgeRet.first);
            boost::put(edge_property_expression_t(), m_Graph, edgeRet.first, 
                double(vertexRec[childVertex]) * expr);
        }        
    }
    return vertex;
}

bool DerivativeGraph::BuildSubgraph(const DervGraphVertex& target, 
                                    const DervGraphVertex& source, 
                                    DervGraph&             subgraph,
                                    std::unordered_map<DervGraphVertex, bool> &vertMap) {
    if(vertMap.find(target) != vertMap.end()) {
        return vertMap[target];
    }
        
    if (source == target) {
        vertMap[target] = true;
        return true;
    }

    bool isSubgraph = false;
    for(auto it = boost::in_edges(target, m_Graph); it.first != it.second; it.first++) {
        auto edge = *it.first;
        if(BuildSubgraph(boost::source(edge, m_Graph), source, subgraph, vertMap)) {
            boost::add_edge(boost::source(edge, m_Graph), 
                            boost::target(edge, m_Graph),
                            boost::get(edge_property_expression_t(), m_Graph, edge),
                            subgraph);
            isSubgraph = true;
        }
    }
    
    vertMap[target] = isSubgraph;
    return isSubgraph;
}

void DerivativeGraph::FactorSubgraphs(std::vector<DerivativeSubgraph>& subgraphs) {
    DervGraphVertex vert0, vert1;
    std::vector<int> subgraphIdList;
    std::vector<bool> isDomList;
    while(GetFactorSubgraph(subgraphs, vert0, vert1, subgraphIdList, isDomList)) {      
        for (int i = 0; i < (int)subgraphIdList.size(); i++) {
            int subgraphId = subgraphIdList[i];
            bool isDom = isDomList[i];
            DerivativeSubgraph &dSubgraph = subgraphs[subgraphId];
            DervGraph &subgraph = dSubgraph.subgraph;
            std::vector<int> &domsArray = isDom ? dSubgraph.postDoms : dSubgraph.doms;        
            std::set<std::pair<int, int> > edges;
            auto factorExpr = FactorExpr(subgraph, vert0, vert1, edges).second;        
            for(auto edge : edges) {            
                DervGraphVertex v = isDom ? boost::vertex(edge.second, subgraph) : boost::vertex(edge.first, subgraph);
                auto targetVert = isDom ? vert1 : vert0;
                bool domTest = v == targetVert || (int)v == domsArray[v];
                if (!domTest) {
                    for(int domId = v; domsArray[domId] != domId; domId = domsArray[domId]) {
                        if(domsArray[domId] == (int)targetVert) {
                            domTest = true;
                            break;
                        }
                    }            
                }
                if(domTest) {                
                    boost::remove_edge(edge.first, edge.second, subgraph);                
                }
            }                
            boost::add_edge(vert0, vert1, factorExpr, subgraph);            
        }
    }
}

void DerivativeGraph::FactorCommonSubproduct(const std::vector<DerivativeSubgraph> &subgraphs,
                                             std::vector<std::shared_ptr<Expression>> &derivatives) {
    typedef std::pair<std::shared_ptr<Expression>, std::shared_ptr<Expression>> EdgePair;
    std::map<EdgePair, int> edgePairMap;
    std::vector<std::list<std::shared_ptr<Expression>>> paths(subgraphs.size());
    std::multimap<std::shared_ptr<Expression>, 
        std::pair<size_t, std::list<std::shared_ptr<Expression>>::iterator>> itMap;
    // For each subgraph, build a path from source to target and store the expression of the edge in a linked list
    // In addition, build a map that maps from the expression to the iterator of the list
    // Note that a single expression can map to multiple iterators because some paths share the same edge
    for(size_t subgraphId = 0; subgraphId < subgraphs.size(); subgraphId++) {        
        const auto& dSubgraph = subgraphs[subgraphId];
        const DervGraph& subgraph = dSubgraph.subgraph;
        auto source = dSubgraph.source;
        auto target = dSubgraph.target;
        std::shared_ptr<Expression> prevExpr;
        while(source != target) {
            auto it = boost::out_edges(source, subgraph);
            auto e = *it.first;
            std::shared_ptr<Expression> expr = 
                boost::get(edge_property_expression_t(), subgraph, e);            
            paths[subgraphId].push_back(expr);
            auto lastIt = std::prev(paths[subgraphId].end());     
            itMap.insert({expr, {subgraphId, lastIt}});
            if (prevExpr.get() != nullptr) {
                EdgePair edgePair(prevExpr, expr);                
                edgePairMap[edgePair] = edgePairMap[edgePair] + 1;                
            }
            prevExpr = expr;
            source = boost::target(e, subgraph); 
        }
    }
    typedef std::pair<EdgePair, int> EdgePairKV;
    std::set<EdgePairKV, std::function<bool(EdgePairKV, EdgePairKV)> > largest(
        [](const EdgePairKV& e0, const EdgePairKV& e1) {
            return e0.second != e1.second ? e0.second < e1.second : e0.first < e1.first;
        });
    for(auto& kv : edgePairMap) {
        largest.insert(kv);
    }    
    while(!largest.empty()) {            
        auto last = std::prev(largest.end());
        auto edgePair = (*last).first;
        auto mergedExpr = edgePair.first * edgePair.second;
        bool found = false;
        for (auto ranges = itMap.equal_range(edgePair.first); ranges.first != ranges.second; ) {     
            auto removeEdgePair = [&](const EdgePair &edgePair) {
                largest.erase(*(edgePairMap.find(edgePair)));                
                edgePairMap[edgePair] = edgePairMap[edgePair] - 1;
                if(edgePairMap[edgePair] > 0) {
                    largest.insert(*(edgePairMap.find(edgePair)));
                } else {
                    edgePairMap.erase(edgePair);
                }                      
            };
            auto insertEdgePair = [&](const EdgePair &edgePair) {
                auto it = edgePairMap.find(edgePair);
                if (it != edgePairMap.end()) {
                    largest.erase(*it);
                }
                edgePairMap[edgePair] = edgePairMap[edgePair] + 1;
                largest.insert(*(edgePairMap.find(edgePair)));            
            };
            auto itMapKV = *ranges.first;
            auto listItPair = itMapKV.second;
            size_t listId = listItPair.first;            
            auto &list = paths[listId];
            auto listIt = listItPair.second;
            auto listItNext = std::next(listIt);
            ranges.first++;
            if (listItNext == list.end()) {                
                continue;
            }
            if (*listIt == edgePair.first && *listItNext == edgePair.second) {                
                found = true;
                // merge listIt & listItNext (update itMap simultaneously)
                // erase listIt & listItNext
                itMap.erase(std::prev(ranges.first));
                for (auto nextRanges = itMap.equal_range(*listItNext); 
                        nextRanges.first != nextRanges.second; nextRanges.first++) {
                    auto itMapKV = *nextRanges.first;
                    auto listItPair = itMapKV.second;
                    if (listItPair.first == listId && listItPair.second == listItNext) {
                        itMap.erase(nextRanges.first);
                        break;
                    }
                }
                auto newIt = list.erase(listIt, std::next(listItNext));
                // insert merged list before listIt
                newIt = list.insert(newIt, mergedExpr);  
                itMap.insert({mergedExpr, {listId, newIt}});
                // This can be faster...?
                ranges = itMap.equal_range(edgePair.first);
                
                // update EdgePair statistics
                removeEdgePair(edgePair);       
                if (newIt != list.begin()) {                    
                    auto newItPrev = std::prev(newIt);
                    EdgePair oldEdgePair(*newItPrev, edgePair.first);
                    removeEdgePair(oldEdgePair);                    
                    EdgePair newEdgePair(*newItPrev, *newIt);
                    insertEdgePair(newEdgePair);                    
                }
                auto newItNext = std::next(newIt);                                
                if (newItNext != list.end()) {
                    EdgePair oldEdgePair(edgePair.second, *newItNext);
                    removeEdgePair(oldEdgePair);
                    EdgePair newEdgePair(*newIt, *newItNext);
                    insertEdgePair(newEdgePair);
                }
            } 
        }

        if (!found) {
            throw std::runtime_error("Error in FactorCommonSubproduct");
        }
    }
    
    // After the factoring, the paths should contain only a single edge, which is the product derivative of the path
    derivatives.resize(paths.size());
    for (size_t i = 0; i < derivatives.size(); i++) {
        derivatives[i] = *(paths[i].begin());
    }
}

} //namespace cdstar