#ifndef TREE_HPP
#define TREE_HPP

#include <vector>
#include <algorithm>
#include <iostream>

#include <flags.hpp>
#include <utils.hpp>

// an element of the octree
class Node
{
    public:
    Node(FL_TYPE ax, FL_TYPE ay, FL_TYPE az, FL_TYPE cur_l, FL_TYPE max_l, int depth = 0);
    size_t numElementsInside();
    bool rayIntersection(
        FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
        FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
        size_t *intersecting_nodes, int &idx,
        bool normalized = false
    );
    
    std::vector<Node *> subnodes;
    std::vector<FL_TYPE> vertices;
    size_t node_id;
    std::vector<size_t> elements;
    bool is_leaf;
    size_t num_all_contained_elements;
    int depth;

    static size_t node_count;
    static std::vector<Node *> all_nodes;
    static std::vector<std::vector<Node *>> node_as_depths;
};


void fillTree(FL_TYPE *element_nodes, size_t num_elements);

void flattenTree(
    std::vector<FL_TYPE> &vec,
    std::vector<FL_TYPE> &vec_pos
);
#endif