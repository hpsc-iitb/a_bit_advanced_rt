#include <tree.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

Node::Node(
    FL_TYPE ax, FL_TYPE ay, FL_TYPE az, FL_TYPE cur_l, FL_TYPE max_l
)
{
    if(cur_l < max_l)
    {
        return;
    }

    if(max_l > cur_l / 2)
    {
        this->is_leaf = true;
    }
    else
    {
        this->is_leaf = false;        /* code */
    }
    
    this->node_id = Node::node_count;
    Node::node_count++;
    /*
     e---f
    /|  /|
    a-g-bh
    |/  |/
    c---d
    */
   // a
   this->vertices.resize(24, 0);

//    this->vertices[0] = ax;
//    this->vertices[1] = ay;
//    this->vertices[2] = az;
    for(size_t _i = 0; _i < this->vertices.size()/3; _i++)
    {
        this->vertices[_i*3 + 0] = ax + (_i%2)*cur_l;
        this->vertices[_i*3 + 1] = ay - cur_l * ((_i%4)>1?1:0);
        this->vertices[_i*3 + 2] = az + cur_l * ((_i/4)>0?1:0);

        // std::cout << this->vertices[_i*3 + 0] << " "\
        //     << this->vertices[_i*3 + 1] << " "\
        //     << this->vertices[_i*3 + 2] << "\n";
        Node *n =  new Node(
            ax + (_i%2)*cur_l,
            ay - cur_l * ((_i%4)>1?1:0),
            az + cur_l * ((_i/4)>0?1:0),
            cur_l / 2,
            max_l
        );
        this->subnodes.push_back(n);
    }
}

size_t Node::elementsInside()
{
    if(this->is_leaf)
    {
        return this->elements.size();
    }

    size_t num_all_elem = 0;
    for(size_t _i = 0; _i < this->subnodes.size(); _i++)
    {
        num_all_elem += this->subnodes[_i]->elements.size();
    }
    return num_all_elem;
}