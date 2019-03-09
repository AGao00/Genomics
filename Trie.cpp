#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <list>
#include <cctype>
#include <iostream>

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
        if (key == "")
            return;
        
            // keeps track of which character in key is being processed
        int index = 0;
        Node* last = pathFound(key, root, index);
        
            // if the path found is exactly the same as key, push value into that Node's vector of ValueTypes
        if (index >= key.size()) {
            last->values.push_back(value);
            return;
        }
        
        while (index < key.size()) {
                // create a new Node with assocated val of first unmatched letter in key, with parent = par
            Node* temp = new Node;
            temp->associated = toupper(key[index]);
            temp->parent = last;

                // push that new Node into par's list of children, set par to that new Node, increment index
            last->children.push_front(temp);
            last = last->children.front();
            index++;
        }
        
            // par now points to the Node that has an associated value of the last letter in key
        last->values.push_back(value);
    }
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const {
        std::vector<ValueType> results;
        
            // if only want exact matches, with no SNiPs
        if (exactMatchOnly) {
                // find last Node that has a match with key
            int index = 0;
            Node* p = pathFound(key, root, index);
            std::vector<ValueType> val = p->values;
            
                // if key is found in Trie, add all values of that Node to results
            if (index >= key.size())
                results.insert(results.end(), val.begin(), val.end());
            
                // return result -> result will be empty if there are no matches
            return results;
        }
        
        // if do want SNiPs
        
        char bases[4] = {'A', 'C', 'T', 'G'};
        for (int i = 1; i < key.size(); i++) {
            for (int j = 0; j < 4; j++) {
                // make SNiPs
                std::string newKey = key.substr(0, i) + bases[j] + key.substr(i+1, key.size());
                
                    // if not a repeat of original key, call find again for newKey with exactMatchOnly as true
                if (newKey != key) {
                    std::vector<ValueType> vals = find(newKey, true);
                    
                        // insert values into overall result
                    results.insert(results.end(), vals.begin(), vals.end());
                }
            }
        }
        
            // process original key
        std::vector<ValueType> vals = find(key, true);
        results.insert(results.end(), vals.begin(), vals.end());
        
        return results;
    }
    
    void printTrie() {
        printing(root, "");
    }

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct Node {
        std::string associated;
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
    void insertNodes(Node* r, const std::string& key, int index) {
        if (index >= key.size())
            return;
        Node* temp = new Node;
        temp->associated = key[index];
        temp->parent = r;
        r->children.push_front(temp);
        insertNodes(r->children.front(), key, index+1);
    }
    Node* pathFound(const std::string& key, Node* p, int& index) const {
        // returns last Node that has a match with key, and index is either end of key or first unmatched character
        // p should never be nullptr, because root is always initialized to point to empty Node
        
            // if p's associated char is equal to current char in key
        if (p->associated == "" || (index < key.size() && p->associated[0] == key[index])) {
                // if there are no more child Nodes, return this Node
            if (p->children.empty()) {
                if (p->associated[0] == key[index])
                    index++;
                return p;
            }
            
                // look to see if there is anymore matches later on
            if (p->associated != "")
                index++;
            for (auto i = p->children.begin(); i != p->children.end(); i++) {
                if ((*i)->associated[0] == key[index])
                    return pathFound(key, (*i), index);
            }
            return p;
        }
        
            // if chars don't match or index exceeds key size, return parent Node
        return p->parent;
    }
    
    void printing(Node* p, std::string tabs) {
        if (p == nullptr) {
            std::cout << std::endl;
            return;
        }
        std::cout << tabs << p->associated << ": ";
        for (int i = 0; i < p->values.size(); i++)
            std::cout << p->values[i] << " ";
        std::cout << std::endl;
        
        for(auto j = p->children.begin(); j != p->children.end(); j++)
            printing(*j, tabs + "  ");
    }
};

#endif // TRIE_INCLUDED
