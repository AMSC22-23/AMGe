template<class T>
class FunctionNode {		//abstract class for mapping
public:
    FunctionNode() = default;
    virtual T getValue(int i,int j)  = 0;
    virtual ~FunctionNode() = default;
};