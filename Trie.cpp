#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <list>

template<typename ValueType>
class Trie
{
public:
    Trie() {
        // empty root Node with no children
        root = new Node;
        root->associated = "";
        root->parent = nullptr;
    }
    ~Trie() { deleteAllNodes(root); }
    void reset() {
        deleteAllNodes(root);
        root = new Node;
        root->associated = "";
        root->parent = nullptr;
    }
    void insert(const std::string& key, const ValueType& value) {
        // keeps track of which character in key is being processed
        int index = 0;
        Node* last = pathFound(key, root, index);
        
            // if the path found is exactly the same as key, push value into that Node's vector of ValueTypes
        if (index >= key.size()) {
            last->values.push_back(value);
            return;
        }
        
        Node* par = last;
        while (index < key.size()) {
                // create a new Node with assocated val of first unmatched letter in key, with parent = par
            Node* temp = new Node;
            temp->associated = key[index];
            temp->parent = par;
            
                // push that new Node into par's list of children, set par to that new Node, increment index
            par->children.push_front(temp);
            par = par->children.front();
            index++;
        }
        
            // par now points to the Node that has an associated value of the last letter in key
        par->values.push_back(value);
    }

    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const {
        std::vector<ValueType> results;
        
        // IMPLEMENT
        
        return results;
    }

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct Node {
        char associated;
        std::vector<ValueType> values;
        std::list<Node*> children;
        Node* parent;
    };
    Node* root;
    void deleteAllNodes(Node* temp) {
        // frees all allocated memory
        if (temp == nullptr)
            return;
        for (auto p = temp->children.begin(); p != temp->children.end(); p++)
            deleteAllNodes(*p);
        delete temp;
    }
    Node* pathFound(const std::string& key, Node* p, int& index) {
        // returns last Node that has a match with key, and index is either end of key or first unmatched character
        // p should never be nullptr, because root is always initialized to point to empty Node
        
            // if p's associated char is equal to current char in key
        if (p->associated == "" || (index < key.size() && p->associated == key[index])) {
            
                // if there are no more child Nodes, return this Node
            if (p->children.empty)
                return p;
            
                // look to see if there is anymore matches later on
            for (auto i = p->children.begin(); i != p->children.end(); i++)
                return pathFound(key, *i, index++);
        }
        
            // if chars don't match or index exceeds key size, return parent Node
        return p->parent;
    }
    
};

#endif // TRIE_INCLUDED
