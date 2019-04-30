#include <tree.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

#include <utils.hpp>

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
    Node::all_nodes.push_back(this);
    /*
     e---f
    /|  /|
    a-g-bh
    |/  |/
    c---d
    */
   // a

//    std::cout << "created cube with {a}: " <<\
//             ax << " "\
//             << ay << " "\
//             << az <<\
            ", length: " << cur_l << "\n";
   this->vertices.resize(24, 0);

//    this->vertices[0] = ax;
//    this->vertices[1] = ay;
//    this->vertices[2] = az;
    for(size_t _i = 0; _i < this->vertices.size()/3; _i++)
    {
        this->vertices[_i*3 + 0] = ax + (_i%2)*cur_l;
        this->vertices[_i*3 + 1] = ay - cur_l * ((_i%4)>1?1:0);
        this->vertices[_i*3 + 2] = az + cur_l * ((_i/4)>0?1:0);

        Node *n =  new Node(
            ax + 0.5*(_i%2)*cur_l,
            ay - 0.5*cur_l * ((_i%4)>1?1:0),
            az + 0.5*cur_l * ((_i/4)>0?1:0),
            cur_l / 2,
            max_l
        );
        this->subnodes.push_back(n);
    }
}

size_t Node::numElementsInside()
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

bool Node::rayIntersection(
    FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
    FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
    size_t *intersecting_nodes, int &idx
)
{
    if(this->numElementsInside() == 0)
    {
        // node doesn't have any element inside, no intersections
        return false;
    }

    // TODO: optimize this normalization
    // normalize();
    normalize(rdx, rdy, rdz, rdx, rdy, rdz);
    rdx = (!rdx)?1e-8:rdx;
    rdy = (!rdy)?1e-8:rdy;
    rdz = (!rdz)?1e-8:rdz;
    
    FL_TYPE swap_tmp;

    FL_TYPE txmin = (this->vertices[0] - rox) / rdx; // vertex a
    FL_TYPE txmax = (this->vertices[3] - rox) / rdx; // vertex b
    
    FL_TYPE tymin = (this->vertices[7] - roy) / rdy; // vertex c
    FL_TYPE tymax = (this->vertices[1] - roy) / rdy; // vertex a

    FL_TYPE tzmin = (this->vertices[2] - roz) / rdz; // vertex a
    FL_TYPE tzmax = (this->vertices[14] - roz) / rdz; // vertex e


    // account for negatives
    if(rdx < 0)
    {
        swap_tmp = txmax;
        txmax = txmin;
        txmin = swap_tmp;
    }
    if(rdy < 0)
    {
        swap_tmp = tymax;
        tymax = tymin;
        tymin = swap_tmp;
    }
    if(rdz < 0)
    {
        swap_tmp = tzmax;
        tzmax = tzmin;
        tzmin = swap_tmp;
    }

    FL_TYPE tmin = std::max(
        txmin, std::max(tymin, tzmin)
    );

    FL_TYPE tmax = std::min(
        txmax, std::min(tymax, tzmax)
    );

    if(tmin <= tmax)
    {
        // there exists a parameter t for which ray intersects block
        if(this->is_leaf)
        {
            // last node, add own id to intersecting_nodes
            intersecting_nodes[idx++] = this->node_id;
            return true;
        }
        else
        {
            bool retval = false;
            for(size_t _k = 0; _k < this->subnodes.size(); _k++)
            {
                // std::cout << "in rayintersection loop, processing: " << _k << "\n";
                retval = retval | this->subnodes[_k]->rayIntersection(
                    rox, roy, roz, rdx, rdy, rdz, intersecting_nodes, idx
                );
            }
            return retval;
        }
    }
    else
    {
        return false;
    }
    
}

size_t Node::node_count = 0;
std::vector<Node *> Node::all_nodes(0);

void fillTree(FL_TYPE *element_nodes, size_t num_elements)
{
    
}