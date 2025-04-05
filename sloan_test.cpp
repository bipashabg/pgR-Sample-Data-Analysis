#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <map>
#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/sloan_ordering.hpp>

typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    boost::property<boost::vertex_color_t, boost::default_color_type,
        boost::property<boost::vertex_degree_t, int,
            boost::property<boost::vertex_priority_t, double> > >,
    boost::property<boost::edge_weight_t, double>
> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertices_size_type size_type;

struct EdgeData {
    int id;
    int source;
    int target;
    double cost;
};

std::pair<Vertex, Vertex> find_pseudoperipheral_pair(const Graph& g) {
    auto vp = vertices(g);
    
    std::vector<Vertex> vertices_vec;
    for (; vp.first != vp.second; ++vp.first) {
        vertices_vec.push_back(*vp.first);
    }
    
    if (vertices_vec.empty()) {
        return {0, 0};
    }
    
    Vertex current = vertices_vec[0];
    
    std::vector<int> distances(num_vertices(g), -1);
    std::vector<Vertex> queue;
    
    distances[current] = 0;
    queue.push_back(current);
    
    size_t front = 0;
    while (front < queue.size()) {
        Vertex v = queue[front++];
        
        for (auto ep = out_edges(v, g); ep.first != ep.second; ++ep.first) {
            Vertex target = boost::target(*ep.first, g);
            if (distances[target] == -1) {
                distances[target] = distances[v] + 1;
                queue.push_back(target);
            }
        }
    }
    
    Vertex end = current;
    int max_dist = 0;
    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] > max_dist) {
            max_dist = distances[i];
            end = i;
        }
    }
    
    std::fill(distances.begin(), distances.end(), -1);
    queue.clear();
    front = 0;
    
    distances[end] = 0;
    queue.push_back(end);
    
    while (front < queue.size()) {
        Vertex v = queue[front++];
        
        for (auto ep = out_edges(v, g); ep.first != ep.second; ++ep.first) {
            Vertex target = boost::target(*ep.first, g);
            if (distances[target] == -1) {
                distances[target] = distances[v] + 1;
                queue.push_back(target);
            }
        }
    }
    
    Vertex start = end;
    max_dist = 0;
    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] > max_dist) {
            max_dist = distances[i];
            start = i;
        }
    }
    
    return {start, end};
}

std::vector<int> apply_sloan_ordering(const Graph& g, const std::map<Vertex, int>& id_map) {
    std::vector<Vertex> sloan_order(num_vertices(g));
    
    auto pseudo_pair = find_pseudoperipheral_pair(g);
    Vertex start_vertex = pseudo_pair.first;
    Vertex end_vertex = pseudo_pair.second;
    
    std::cout << "Using pseudoperipheral vertices for Sloan algorithm:" << std::endl;
    std::cout << "Start vertex index: " << start_vertex << " (ID: " << id_map.at(start_vertex) << ")" << std::endl;
    std::cout << "End vertex index: " << end_vertex << " (ID: " << id_map.at(end_vertex) << ")" << std::endl;
    
    std::vector<int> degree_vec(num_vertices(g));
    auto degree_map = boost::make_iterator_property_map(degree_vec.begin(), get(boost::vertex_index, g));
    
    auto vp = vertices(g);
    for (; vp.first != vp.second; ++vp.first) {
        degree_map[*vp.first] = out_degree(*vp.first, g);
    }
    
    std::vector<double> priority_vec(num_vertices(g));
    auto priority_map = boost::make_iterator_property_map(priority_vec.begin(), get(boost::vertex_index, g));
    
    std::vector<boost::default_color_type> color_vec(num_vertices(g), boost::white_color);
    auto color_map = boost::make_iterator_property_map(color_vec.begin(), get(boost::vertex_index, g));
    
    sloan_ordering(
        g,
        std::back_inserter(sloan_order),
        color_map,
        degree_map,
        priority_map,
        start_vertex,
        end_vertex
    );
    
    std::vector<int> result;
    for (const auto& v : sloan_order) {
        result.push_back(id_map.at(v));
    }
    
    return result;
}

int calculate_bandwidth(const std::vector<EdgeData>& edges, const std::vector<int>& order) {
    std::map<int, int> position_map;
    for (size_t i = 0; i < order.size(); ++i) {
        position_map[order[i]] = i;
    }
    
    int max_bandwidth = 0;
    for (const auto& edge : edges) {
        auto src_pos = position_map.find(edge.source);
        auto tgt_pos = position_map.find(edge.target);
        
        if (src_pos != position_map.end() && tgt_pos != position_map.end()) {
            int bandwidth = std::abs(src_pos->second - tgt_pos->second);
            max_bandwidth = std::max(max_bandwidth, bandwidth);
        }
    }
    
    return max_bandwidth;
}

int calculate_profile(const std::vector<EdgeData>& edges, const std::vector<int>& order) {
    std::map<int, int> position_map;
    for (size_t i = 0; i < order.size(); ++i) {
        position_map[order[i]] = i;
    }
    
    int total_profile = 0;
    for (const auto& edge : edges) {
        auto src_pos = position_map.find(edge.source);
        auto tgt_pos = position_map.find(edge.target);
        
        if (src_pos != position_map.end() && tgt_pos != position_map.end()) {
            total_profile += std::abs(src_pos->second - tgt_pos->second);
        }
    }
    
    return total_profile;
}

int main() {
    std::vector<EdgeData> edges = {
        {1, 1, 2, 1.0},
        {2, 2, 3, 1.0},
        {3, 3, 4, 1.0},
        {4, 4, 5, 1.0},
        {5, 1, 6, 1.0},
        {6, 6, 7, 1.0},
        {7, 7, 8, 1.0},
        {8, 8, 9, 1.0},
        {9, 9, 17, 1.0},
        {10, 2, 7, 1.0},
        {11, 3, 8, 1.0},
        {12, 4, 9, 1.0},
        {13, 5, 17, 1.0},
        {14, 2, 1, 1.0},
        {15, 3, 2, 1.0},
        {16, 4, 3, 1.0},
        {17, 5, 4, 1.0},
        {18, 6, 1, 1.0},
        {19, 7, 6, 1.0},
        {20, 8, 7, 1.0},
        {21, 9, 8, 1.0},
        {22, 17, 9, 1.0},
        {23, 7, 2, 1.0},
        {24, 8, 3, 1.0},
        {25, 9, 4, 1.0},
        {26, 17, 5, 1.0}
    };
    
    Graph g;
    
    std::map<int, Vertex> vertex_map;
    std::map<Vertex, int> id_map;
    
    std::set<int> all_vertices;
    for (const auto& e : edges) {
        all_vertices.insert(e.source);
        all_vertices.insert(e.target);
    }
    
    for (const auto& v : all_vertices) {
        Vertex vd = add_vertex(g);
        vertex_map[v] = vd;
        id_map[vd] = v;
    }
    
    for (const auto& e : edges) {
        Vertex src = vertex_map[e.source];
        Vertex tgt = vertex_map[e.target];
        
        add_edge(src, tgt, e.cost, g);
    }
    
    std::cout << "pgRouting Sample Data - Sloan Ordering Test" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Graph has " << all_vertices.size() << " vertices and " 
              << num_edges(g) << " edges." << std::endl;
    
    std::vector<int> natural_order(all_vertices.begin(), all_vertices.end());
    
    std::vector<int> sloan_order = apply_sloan_ordering(g, id_map);
    
    int natural_bandwidth = calculate_bandwidth(edges, natural_order);
    int sloan_bandwidth = calculate_bandwidth(edges, sloan_order);
    
    int natural_profile = calculate_profile(edges, natural_order);
    int sloan_profile = calculate_profile(edges, sloan_order);
    
    std::cout << "\nNatural ordering bandwidth: " << natural_bandwidth << std::endl;
    std::cout << "Natural ordering profile: " << natural_profile << std::endl;
    
    std::cout << "\nSloan ordering bandwidth: " << sloan_bandwidth << std::endl;
    std::cout << "Sloan ordering profile: " << sloan_profile << std::endl;
    
    if (sloan_bandwidth < natural_bandwidth) {
        std::cout << "\nBandwidth reduction: " 
                  << (100.0 * (natural_bandwidth - sloan_bandwidth) / natural_bandwidth) 
                  << "%" << std::endl;
    } else {
        std::cout << "\nWarning: Bandwidth increased by " 
                  << (100.0 * (sloan_bandwidth - natural_bandwidth) / natural_bandwidth) 
                  << "%" << std::endl;
    }
    
    if (sloan_profile < natural_profile) {
        std::cout << "Profile reduction: " 
                  << (100.0 * (natural_profile - sloan_profile) / natural_profile) 
                  << "%" << std::endl;
    } else {
        std::cout << "Warning: Profile increased by " 
                  << (100.0 * (sloan_profile - natural_profile) / natural_profile) 
                  << "%" << std::endl;
    }
    
    return 0;
}
