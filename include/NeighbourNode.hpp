#ifndef __NEIGHBOURNODE_HPP__
#define __NEIGHBOURNODE_HPP__

struct NeighbourNode {
    public:
    NeighbourNode(){}
    NeighbourNode(int index_val, double weigh, bool cardinal_bool){
        index = index_val;
        weight = weight;
        cardinal = cardinal_bool;
    }
    
    int index;
    double weight;
    int cardinal;
};

#endif