#ifndef __TRIPLET_HPP__
#define __TRIPLET_HPP__

struct Node {
    public:
    Node(){}
    Node(int index_val, double x_val, double y_val, char c){
        index_node = index_val;
        x = x_val;
        y = y_val;
        type = c;
    }
    
    int index_node;
    double x;
    double y;
    int num_neighbours;
    char type;
};

#endif