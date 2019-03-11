#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <cctype>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
    void print();
private:
    string m_name;
    string m_sequence;
    unsigned int m_length;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
    m_length = static_cast<int>(sequence.length());
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    std::string line;
    std::string name;
    std::string genome;
    
    while (getline(genomeSource, line)) {
        if (line[0] == '>') {
                // if start of new sequence, create Genome object with previous sequence
            if (name != "" && genome != "")
                genomes.push_back(Genome(name, genome));
            
                // if invalid state, return false
            if (line.size() == 1 || !isalnum(line[1]))
                return false;
            
                // set name to new name, and genome to empty
            name = line.substr(1, line.length()-1);
            genome = "";
            continue;
        }
            // if no name for sequence
        if (name == "")
            return false;
        
            // check each character to see if they are ACTGN
        for (int i = 0; i < line.size(); i ++) {
            char c = toupper(line[i]);
            line[i] = c;
            switch (c) {
                case 'A':
                case 'G':
                case 'T':
                case 'C':
                case 'N':
                    break;
                default:
                    return false;
            }
        }
            // append line to genome
        genome += line;
    }
    
    if (name == "" || genome == "")
        return false;
    
        // push last Genome object into genomes
    genomes.push_back(Genome(name, genome));
    
    return true;
}

int GenomeImpl::length() const
{
    return m_length;
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
        // if length causes out of bounds, or position/length are invalid values, return false
    if (position + length > m_length || position < 0 || length < 0)
        return false;
    
        // set fragment to substring and return true
    fragment = m_sequence.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
