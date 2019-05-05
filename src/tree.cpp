#include <tree.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

#include <utils.hpp>

Node::Node(
    FL_TYPE ax, FL_TYPE ay, FL_TYPE az, FL_TYPE cur_l, FL_TYPE max_l, int depth
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
        this->is_leaf = false;
    }
    this->node_id = Node::node_count;
    Node::node_count++;
    Node::all_nodes.push_back(this);
    Node::node_as_depths.at(depth).push_back(this);
    this->depth = depth;
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
   this->num_all_contained_elements = 0;

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
            max_l, this->depth + 1
        );
        this->subnodes.push_back(n);
    }
}

size_t Node::numElementsInside()
{
    if(this->is_leaf)
    {
        this->num_all_contained_elements = this->elements.size();
    }
    else
    {
        size_t num_all_elem = 0;
        for(size_t _i = 0; _i < this->subnodes.size(); _i++)
        {
            num_all_elem += this->subnodes[_i]->numElementsInside();
        }
        this->num_all_contained_elements = num_all_elem;
    }
    return this->num_all_contained_elements;
}

bool Node::rayIntersection(
    FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
    FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
    size_t *intersecting_nodes, int &idx,
    bool normalized
)
{
    if(this->num_all_contained_elements == 0)
    {
        // node doesn't have any element inside, no intersections
        return false;
    }

    // TODO: optimize this normalization
    // normalize();
    rdx = (!rdx)?1e-8:rdx;
    rdy = (!rdy)?1e-8:rdy;
    rdz = (!rdz)?1e-8:rdz;
    if(!normalized)
    {
        normalize(rdx, rdy, rdz, rdx, rdy, rdz);
    }
    
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
                retval = retval | this->subnodes[_k]->rayIntersection(
                    rox, roy, roz, rdx, rdy, rdz, intersecting_nodes, idx,
                    true
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
std::vector<std::vector <Node *>> Node::node_as_depths(max_depth + 1);

void fillTree(FL_TYPE *element_nodes, size_t num_elements)
{
    FL_TYPE xmax, xmin, ymax, ymin, zmax, zmin;
    size_t eid;
    for (size_t _i = 0; _i < Node::node_count; _i++)
    {
        Node *n = Node::all_nodes.at(_i);
        // don't store if not leaf
        if(n->is_leaf)
        {
            xmax = n->vertices[3]; // b
            xmin = n->vertices[0]; // a
            ymax = n->vertices[1]; // a
            ymin = n->vertices[7]; // c
            zmax = n->vertices[14]; // e
            zmin = n->vertices[2]; // e
            for (size_t _j = 0; _j < num_elements; _j++)
            {
                eid = _j * element_size;
                if(
                    (element_nodes[eid] <= xmax && element_nodes[eid] >= xmin\
                        && element_nodes[eid+1] <= ymax && element_nodes[eid+1] >= ymin\
                        && element_nodes[eid+2] <= zmax && element_nodes[eid+2] >= zmin) ||
                        (element_nodes[eid+3] <= xmax && element_nodes[eid+3] >= xmin\
                        && element_nodes[eid+4] <= ymax && element_nodes[eid+4] >= ymin\
                        && element_nodes[eid+5] <= zmax && element_nodes[eid+5] >= zmin) ||
                        (element_nodes[eid+6] <= xmax && element_nodes[eid+6] >= xmin\
                        && element_nodes[eid+7] <= ymax && element_nodes[eid+7] >= ymin\
                        && element_nodes[eid+8] <= zmax && element_nodes[eid+8] >= zmin)
                )
                {
                    n->elements.push_back(_j);
                }
            }
            
        }
    }
}

void flattenTree(
    std::vector<FL_TYPE> &vec,
    std::vector<FL_TYPE> &vec_pos
    )
{
    vec.resize(0);
    vec_pos.resize(Node::all_nodes.size(), -1);
    /*
    id: 1, vertices: 24, is_leaf: 1, num_contained_elems: 1,
    (if not a leaf){ids of children: 8},
    (if leaf){idx of contained elements}
    */
    for(size_t _i = 0; _i < Node::node_as_depths.size(); _i++)
    {
        // std::cout << Node::node_as_depths.at(_i).size() << "\n";
        for(size_t _j = 0; _j < Node::node_as_depths.at(_i).size(); _j++)
        {
            Node *n = Node::node_as_depths.at(_i).at(_j);
            vec_pos[n->node_id] = vec.size();
            
            //push id
            vec.push_back(n->node_id);
            for(size_t _k = 0; _k < 24; _k++)
            {
                vec.push_back(n->vertices[_k]);
            }
            vec.push_back(n->is_leaf);
            vec.push_back(n->num_all_contained_elements);
            if(!n->is_leaf)
            {
                for(size_t _k = 0; _k < n->subnodes.size(); _k++)
                {
                    vec.push_back(n->subnodes.at(_k)->node_id);
                }                
            }
            else
            {
                for(size_t _k = 0; _k < n->elements.size(); _k++)
                {
                    vec.push_back(n->elements.at(_k));
                }
            }
        }
    }
}