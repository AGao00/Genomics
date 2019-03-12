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
using intPair = pair<int, int>;
using intPairList = list<intPair>;

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
    
    vector<DNAMatch> findGenomesHelper(const string& fragment, int minimumLength) const;
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
        p.first = genome.name();
        p.second = index;
        library->insert(sequence, p);
        index++;
    }
    
    if (index < genome.length() && genome.extract(index, genome.length()-index, sequence)) {
        Pair p;
        p.first = genome.name();
        p.second = index;
        library->insert(sequence, p);
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const {
    if (fragment.size() < minimumLength || minimumLength < minimumSearchLength())
        return false;
    
        // stores all results on matches
    vector<DNAMatch> results;
    
        // if exactMatchOnly
    if (exactMatchOnly) {
        auto x = findGenomesHelper(fragment, minimumLength);
        results.insert(results.end(), x.begin(), x.end());
    }
    else {
            // calculate all possible SNiPs
        char bases[4] = { 'A', 'C', 'G', 'T' };
        for (int i = 1; i < fragment.size(); i++) {
            for (int j = 0; j < 4; j++) {
                auto x = findGenomesHelper(fragment.substr(0, i) + bases[j] + fragment.substr(i+1), minimumLength);
                results.insert(results.end(), x.begin(), x.end());
            }
        }
    }
    
        // if there are no matches
    if (results.empty())
        return false;
    
        // set matches to results
    matches = results;
    
    return true;  // This compiles, but may not be correct
}

vector<DNAMatch> GenomeMatcherImpl::findGenomesHelper(const string& fragment, int minimumLength) const {
    // preconditions: fragment and minimumLength are appropriate lengths, called only when exactMatchOnly is true
    
    vector<DNAMatch> results;
    
    int length = minimumSearchLength();
        // vector of pairs of ints hold potential matches, and the pairs hold starting and ending position of match
    unordered_map<string, intPairList> genomes;
    vector<Pair> res = library->find(fragment.substr(0, length), true);
    
        // put all Pairs in res into genomes, so can keep track of number of possible matches
    for (int i = 0; i < res.size(); i++) {
        string name = res[i].first;
        int pos = res[i].second;
        bool repeat = false;
        if (genomes.find(name) != genomes.end()) {
            auto list = genomes[name];
            for (auto pointer = list.begin(); pointer != list.end(); pointer++) {
                if ((*pointer).second > pos)
                    repeat = true;
            }
        }
        if (!repeat)
            genomes[name].push_back(intPair(pos, pos+length));
    }
    
    for (int i = 1; i < fragment.size(); i++) {
            // make sure length won't exceed bounds of string
        if (i+length > fragment.size())
            length = static_cast<int>(fragment.size())-i;
        
            // find new subsequence of fragment
        res = library->find(fragment.substr(i, length), true);
        
            // go through all pairs found in library
        for (int j = 0; j < res.size(); j++) {
            int pos = res[j].second;
            string name = res[j].first;
            
                // if name not in original list of possible matches, discard as viable match
            if (genomes.find(name) == genomes.end())
                continue;
            
            auto v = genomes[name];
            for (auto k = v.begin(); k != v.end();) {
                auto p = *k;
                
                    // if found match of subsequence links up to original match, extend position of last character matching
                if (p.first+i == pos) {
                    p.second = pos+length;
                    k++;
                }
                    // if potential match's length is less than minimumLength and will never be enlongated, delete
                else if (p.second - p.first < minimumLength)
                    k = v.erase(k);
            }
        }
    }
    
        // genomes should be filled with Genome names and the starting/ending position of the matches
    for (auto match : genomes) {
            // for all matches found in that one genome, get starting position and length
        for (auto i = match.second.begin(); i != match.second.end(); i++) {
            DNAMatch m;
            m.genomeName = match.first;
            m.position = (*i).first;
            m.length = (*i).second - m.position;
            if (m.length < minimumLength)
                continue;
            results.push_back(m);
        }
    }
    
    return results;
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
