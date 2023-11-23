#ifndef __FUNCTION_NODE_HPP__
#define __FUNCTION_NODE_HPP__


template<class T>
class FunctionNode {		//abstract class for mapping
public:
    FunctionNode()                   = default;
    virtual ~FunctionNode()          = default;
    virtual T getValue(int i, int j) = 0;
};


#endif
