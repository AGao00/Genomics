#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <unordered_map>
using namespace std;

using Pair = pair<string, int>;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    Trie<Pair>* library;
    int m_minLength;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minLength = minSearchLength;
    library = new Trie<Pair>;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
        // index will keep track of where the next subsequence is, and sequence contains current subsequence
    int index = 0;
    string sequence;
    
    while (genome.extract(index, m_minLength, sequence)) {
        Pair p;
        p.first = sequence;
        p.second = index;
        library->insert(sequence, p);
        index++;
    }
    
    if (index < genome.length() && genome.extract(index, genome.length()-index, sequence)) {
        Pair p;
        p.first = sequence;
        p.second = index;
        library->insert(sequence, p);
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if (fragment.size() < minimumLength || minimumLength < minimumSearchLength())
        return false;
    
    vector<DNAMatch> results;
    
    int index = 0, length = minimumSearchLength(), remaining_length = static_cast<int>(fragment.size())-length;
        // res is vector of all pairs that have path in Trie library equal to first subsequence of fragment
    vector<Pair> res = library->find(fragment.substr(index, length), exactMatchOnly);
    
        // unordered_map of Genome name mapped to a vector of positions
    unordered_map<string, vector<int>> possible_matches;
    for (int i = 0; i < res.size(); i++) {
        Pair p = res[i];
        possible_matches[p.first].push_back(p.second);
    }
    
        // while there is remaining substrings of fragment not yet processed
    while (remaining_length > 0) {
        index += length;
        if (remaining_length < length)
            length = remaining_length;
        
            // iterate through all found Pairs of new subsequence
        res = library->find(fragment.substr(index, length), exactMatchOnly);
        for (int i = 0; i < res.size(); i++) {
            Pair p = res[i];
            auto vals = possible_matches.find(p.first);
            
                // if Genome name has been found before, add position to vector
            if (vals != possible_matches.end())
                possible_matches[p.first].push_back(p.second);
        }
        
        remaining_length -= length;
    }
    
    // now we have a hashtable of Genome names and positions within those names that are possibly matches of fragment
    for (auto x : possible_matches) {
        string gen = x.first;
        vector<int> values = x.second;
        
    }
    
    
    return false;  // This compiles, but may not be correct
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return false;  // This compiles, but may not be correct
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
