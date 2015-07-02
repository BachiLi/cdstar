#ifndef CDSTAR_DERIVATIVE_GRAPH_H__
#define CDSTAR_DERIVATIVE_GRAPH_H__

#include "expression.h"

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace cdstar {

struct edge_property_expression_t {
    typedef boost::edge_property_tag kind;
};

struct vertex_property_expression_t {
    typedef boost::vertex_property_tag kind;
};

typedef boost::property<vertex_property_expression_t, std::shared_ptr<Expression>> VertexProperty;
typedef boost::property<edge_property_expression_t, std::shared_ptr<Expression>, 
                        boost::property<boost::edge_index_t, int>> EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, 
                              VertexProperty, EdgeProperty> DervGraph;
typedef boost::graph_traits<DervGraph>::vertex_descriptor DervGraphVertex;
typedef boost::graph_traits<DervGraph>::edge_descriptor DervGraphEdge;

struct DerivativeSubgraph {
    DervGraph subgraph;
    DervGraphVertex source;
    DervGraphVertex target;
    std::vector<int> doms;
    std::vector<int> postDoms;
};

class DerivativeGraph {
public:
    DerivativeGraph(const std::vector<ExprPtrPair> &dervExprs);
    std::vector<std::shared_ptr<Expression>> Derivatives() const {
        return m_Derivatives;
    }
private:        
    DervGraphVertex BuildGraph(const std::shared_ptr<Expression> node);
    bool BuildSubgraph(const DervGraphVertex &target, 
                       const DervGraphVertex &source, 
                       DervGraph &subgraph,
                       std::unordered_map<DervGraphVertex, bool> &vertMap);
    bool GetFactorSubgraph(std::vector<DerivativeSubgraph> &subgraphs,                     
                           DervGraphVertex &vert0, 
                           DervGraphVertex &vert1,
                           std::vector<int> &subgraphId,
                           std::vector<bool> &isDom);
    void FactorSubgraphs(std::vector<DerivativeSubgraph> &subgraphs);
    void FactorCommonSubproduct(const std::vector<DerivativeSubgraph> &subgraphs,
                                std::vector<std::shared_ptr<Expression>> &derivatives);
    std::unordered_map<const Expression*, DervGraphVertex> m_Vertices;
        
    DervGraph m_Graph;
    std::vector<std::shared_ptr<Expression>> m_Derivatives;
};

} // namespace cdstar

#endif //CDSTAR_DERIVATIVE_GRAPH_H__