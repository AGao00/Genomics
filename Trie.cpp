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
        // inserts value at last node of path that spells out key
        
        if (key == "")
            return;
        
            // keeps track of which character in key is being processed
        int index = 0;
        Node* last = pathFound(key, root, index); // O(LC) L is key.length(), C is avg # of children
        
            // if the path found is exactly the same as key, push value into that Node's vector of ValueTypes
        if (index >= key.size()) {
            last->values.push_back(value);
            return;
        }
        
        while (index < key.size()) { // O(L) where L is remaining characters in key not found earlier
                // create a new Node with assocated val of first unmatched letter in key, with parent = par
            Node* temp = new Node;
            temp->associated = key[index];
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
        return findSNiPs(key, root, exactMatchOnly);
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
    Node* pathFound(const std::string& key, Node* p, int& index) const {
        // returns last Node that has a match with key, and index is either end of key or first unmatched character
        // p should never be nullptr, because root is always initialized to point to empty Node
        // O(LC) where L is length of key and C is average number of children per Node
        
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
    std::vector<ValueType> findSNiPs(const std::string& key, Node* p, bool exactMatchOnly) const {
            // if do want SNiPs, check all children of first letter
        
        std::vector<ValueType> results;
        
            // if exactMatchOnly, run pathFound to find last Node that matches
        if (exactMatchOnly) {
            int count = 0;
            Node* l = pathFound(key, p, count);
            
                // if entire key matched, add values at that Node
            if (count >= key.size())
                results.insert(results.end(), l->values.begin(), l->values.end());
            return results;
        }
        
            // find SNiPs by starting with first letter that must stay constant
        int count = 0;
        Node* first_letter = pathFound(key.substr(0,1), p, count);
        
            // if key is one letter, no SNiPs
        if (count == key.size()) {
            return first_letter->values;
        }
        
            // store all children of first letter
        std::list<Node*> first_children = first_letter->children;
        
        for (auto it = first_children.begin(); it != first_children.end(); it++) {
                // n is starting Node to be checked, c is similar to counter
            Node* n = *it;
            int counter = count;
            int c = 0;
                // new key to be processed - replaces one character with associated value of current node
            std::string newKey = n->associated[0] + key.substr(counter+1);
            std::vector<ValueType> res;
            
                // if no mismatch, keep iterating through children to find one
            if (newKey == key.substr(counter)) {
                    // if no more children and newkey is one character, add values to results
                if (n->children.empty() && newKey.size() == 1)
                    res = n->values;
                
                    // continue with iterating through key for mismatches
                for (auto it2 = n->children.begin(); it2 != n->children.end(); it2++) {
                    res = findSNiPs(newKey, n, exactMatchOnly);
                }
            }
                // there is a mismatch
            else {
                Node* last = pathFound(newKey, n, c);
                    // if there are no more matches other than mismatch, then not SNiP
                if (last == n && newKey.size() > 1)
                    continue;
                
                counter += c;
                    // if SNiP, add values
                if (counter >= key.size())
                    res = last->values;
            }
                // insert all values found into results
            results.insert(results.end(), res.begin(), res.end());
        }
        return results;
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
