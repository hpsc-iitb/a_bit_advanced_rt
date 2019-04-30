#ifndef TREE_HPP
#define TREE_HPP

#include <vector>
#include <algorithm>
#include <iostream>

#include <flags.hpp>

// an element of the octree
class Node
{
    public:
    Node(FL_TYPE ax, FL_TYPE ay, FL_TYPE az, FL_TYPE cur_l, FL_TYPE max_l);
    size_t elementsInside();
    
    std::vector<Node *> subnodes;
    std::vector<FL_TYPE> vertices;
    size_t node_id;
    std::vector<size_t> elements;
    bool is_leaf;

    static size_t node_count;
};

#endif