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
using PairVector = vector<Pair>;

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
    
    static bool PairComp(Pair x, Pair y);
    bool isAMatch(string& sequence, const string& fragment, int minLength) const;
    bool isASNiP(string& sequence, const string& fragment, int minLength) const;
    unordered_map<string, PairVector> findGenomesHelper(const string& fragment, int minLength) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minLength = minSearchLength;
    m_library = new Trie<Pair>;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    Genome* insert = new Genome(genome);
    genomeLibrary[genome.name()] = insert;
    
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
    
    vector<DNAMatch> m;
        // string is genome name, vector of <string subsequence, int pos>
    unordered_map<string, PairVector> results = findGenomesHelper(fragment, minimumLength);
    
        // find SNiPs
    if (!exactMatchOnly) {
        char bases[5] = { 'A', 'C', 'T', 'G', 'N' };
        for (int i = 1; i < fragment.size(); i++) {
            for (int j = 0; j < 5; j++) {
                string frag = fragment.substr(0, i) + bases[j] + fragment.substr(i+1);
                if (frag == fragment)
                    continue;
                unordered_map<string, PairVector> r = findGenomesHelper(frag, minimumLength);
                
                    // append resulting map into overarching map
                for (auto it = r.begin(); it != r.end(); it++) {
                    pair<string, PairVector> item = *it;
                    results[item.first].insert(results[item.first].end(), item.second.begin(), item.second.end());
                }
            }
        }
    }
    
    for (auto it = results.begin(); it != results.end(); it++) {
        pair<string, PairVector> item = *it;
        if (item.second.empty())
            continue;
        
            // find longest match by heapifying all Paris and getting front()
        make_heap(item.second.begin(), item.second.end(), PairComp);
        Pair x = item.second.front();
        
            // create new DNAMatch and add to m, which will become matches
        DNAMatch toAdd;
        toAdd.genomeName = item.first;
        toAdd.length = x.first.size();
        if (toAdd.length == 0)
            continue;
        toAdd.position = x.second;
        m.push_back(toAdd);
    }
    
    if (m.empty())
        return false;
    
    matches = m;
    
    return true;  // This compiles, but may not be correct
}

unordered_map<string, PairVector> GenomeMatcherImpl::findGenomesHelper(const string& fragment, int minLength) const {
    unordered_map<string, PairVector> matches;
    vector<Pair> potentials = m_library->find(fragment.substr(0, minimumSearchLength()), true);
    
        // find Genomes with matches of first subsequence of fragment
    for (int i = 0; i < potentials.size(); i++) {
        string name = potentials[i].first, extract;
        int pos = potentials[i].second;
        if (genomeLibrary.find(name) == genomeLibrary.end())
            continue;
        if (genomeLibrary.at(name)->extract(pos, static_cast<int>(fragment.size()), extract)) {
                // if rest of extract is a match to fragment, add it to map
            if (isAMatch(extract, fragment, minLength)) {
                Pair toAdd;
                toAdd.first = extract;
                toAdd.second = pos;
                
                matches[name].push_back(toAdd);
            }
        }
    }
    
        // set vector to contain only largest match
    for (auto it = matches.begin(); it != matches.end(); it++) {
        pair<string, PairVector> item = *it;
        make_heap(item.second.begin(), item.second.end(), PairComp);
        Pair longest = item.second.front();
        item.second.clear();
        item.second.push_back(longest);
    }
    
    return matches;
}

bool GenomeMatcherImpl::isAMatch(string& sequence, const string& fragment, int minLength) const {
    // return true if sequence is longer than minLength and ALSO matches with fragment, modifying sequence appropriately
    
        // if sequence is equal to fragment, then true
    if (sequence == fragment)
        return true;
    
    int length = static_cast<int>(fragment.size())-1;
    for (int i = length; i >= minLength; i--) {
            // if some subsequence of sequence is equal to the same subsequence of fragment, then true
        if (sequence.substr(0, i) == fragment.substr(0, i)) {
            sequence = sequence.substr(0, i);
            return true;
        }
    }
    
    return false;
}

bool GenomeMatcherImpl::PairComp(Pair x, Pair y) {
    // if true, x should be before/less than y
    
        // if length of x's extract is longer than length of y
    if (x.first.size() > y.first.size())
        return false;
    
        // if lengths are equal, and x's starting position is earlier than y's
    if (x.first.size() == y.first.size() && x.second < y.second)
        return false;
    
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const {
    
    if (query.length() < fragmentMatchLength || fragmentMatchLength < minimumSearchLength())
        return false;
    
        // string for Genome name, int for number of S matches
    unordered_map<string, int> matches;
    int S = query.length()/fragmentMatchLength;
    
        // for all sequences of length fragmentMatchLength
    for (int i = 0; i < S; i++) {
        vector<DNAMatch> m;
        string fragment;
        query.extract(i*fragmentMatchLength, fragmentMatchLength, fragment);
        
            // find genome names that have matches
        if (!findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, m))
            continue;
        
            // if match, increase count in matches
        for (int i = 0; i < m.size();i++) {
            matches[m[i].genomeName] ++;
        }
    }
    
    vector<GenomeMatch> percentages;
    
        // for all matched genomes, calculate percentage
    for (auto it = matches.begin(); it != matches.end(); it++) {
        GenomeMatch m;
        m.genomeName = (*it).first;
        m.percentMatch = (*it).second/static_cast<double>(S) * 100;
        
            // if percentage is above threshold, push into percentages
        if (m.percentMatch >= matchPercentThreshold)
            percentages.push_back(m);
    }
    
    if (percentages.empty())
        return false;
    
    results = percentages;
    
    return true;  // This compiles, but may not be correct
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
