#ifndef __TRIPLET_HPP__
#define __TRIPLET_HPP__

struct Node {
    public:
    Node(int index_val, double x_val, double y_val){
        index_node = index_val;
        x = x_val;
        y = y_val;
    }
    
    int index_node;
    double x;
    double y;
    int node_neighbours;
};

#endif