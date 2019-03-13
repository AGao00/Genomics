#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <unordered_map>
#include <algorithm>
using namespace std;

using Pair = pair<string, int>;
using intPair = pair<int, int>;
using intPairVector = vector<intPair>;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    Trie<Pair>* m_library;
    int m_minLength;
    unordered_map<string, const Genome*> genomeLibrary;
    
    static bool intPairComp(intPair x, intPair y);
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minLength = minSearchLength;
    m_library = new Trie<Pair>;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    genomeLibrary[genome.name()] = &genome;
    
        // index will keep track of where the next subsequence is, and sequence contains current subsequence
    int index = 0;
    string sequence;
    
    while (genome.extract(index, m_minLength, sequence)) {
        Pair p;
        p.first = genome.name();
        p.second = index;
        m_library->insert(sequence, p);
        index++;
    }
    
    if (index < genome.length() && genome.extract(index, genome.length()-index, sequence)) {
        Pair p;
        p.first = genome.name();
        p.second = index;
        m_library->insert(sequence, p);
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const {
    if (fragment.size() < minimumLength || minimumLength < minimumSearchLength())
        return false;
    
    vector<DNAMatch> results;
    
    if (exactMatchOnly) {
        unordered_map<string, intPairVector> stores;
        int length = minimumSearchLength();
        vector<Pair> res = m_library->find(fragment.substr(0, length), exactMatchOnly);
        
        for (int i = 0; i < res.size(); i++) {
            intPair temp;
            temp.first = res[i].second;
            temp.second = length;
            stores[res[i].first].push_back(temp);
        }
        
        for (int index = 1; index < fragment.size()-length; index++) {
            res = m_library->find(fragment.substr(index, length), exactMatchOnly);
            for (int i = 0; i < res.size(); i++) {
                if (stores.find(res[i].first) == stores.end())
                    continue;
                intPairVector temp = stores.at(res[i].first);
                int pos = res[i].second;
                bool continuation = false;
                for (auto it = temp.begin(); it != temp.end(); it++) {
                    int startpos = (*it).first, seq_len = (*it).second;
                    if (startpos+seq_len+1 >= pos+length) {
                        (*it).second++;
                        continuation = true;
                    }
                }
                if (!continuation) {
                    intPair t;
                    t.first = pos;
                    t.second = length;
                    temp.push_back(t);
                }
            }
        }
        
        for (auto it = stores.begin(); it != stores.end(); it++) {
            vector<intPair> temp = (*it).second;
            make_heap(temp.begin(), temp.end(), GenomeMatcherImpl::intPairComp);
            intPair t = temp.front();
            if (t.second >= minimumLength) {
                DNAMatch m;
                m.genomeName = (*it).first;
                m.length = t.second;
                m.position = t.first;
                results.push_back(m);
            }
        }
    }
    else {
        char bases[5] = { 'A', 'C', 'T', 'G', 'N'};
        vector<DNAMatch> temp;
        for (int i = 1; i < fragment.size(); i++) {
            for (int j = 0; j < 5; j++) {
                string newKey = fragment.substr(0, i) + bases[j] + fragment.substr(i+1);
                if (newKey == fragment)
                    continue;
                if (findGenomesWithThisDNA(newKey, minimumLength, true, temp))
                    results.insert(results.end(), temp.begin(), temp.end());
            }
        }
        if (findGenomesWithThisDNA(fragment, minimumLength, true, temp))
            results.insert(results.end(), temp.begin(), temp.end());
    }
    
    if (results.empty())
        return false;
    
    matches = results;
    
    return true;  // This compiles, but may not be correct
}

bool GenomeMatcherImpl::intPairComp(intPair x, intPair y) {
    // if true, x should be before y
    
        // if length of x is longer than length of y
    if (x.second > y.second)
        return true;
    
        // if lengths are equal, and x's starting position is earlier than y's
    if (x.second == y.second && x.first < y.first)
        return true;
    
    return false;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const {
    
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
