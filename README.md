# Project_4
UCLA CS 32 Project 4 - Genomics

THREE CLASSES:

Trie class
Purpose: Maintain a library of genomes; user can add new genomes
Functionality:
  Create a multimap class with a set key value of string
  When user searches for key, return values associated with key and its SNiPs
  Use string, list classes

Genome class
Purpose:
  User can search library for specific DNA sequence
  Return genomes that contain the DNA sequence
  Return genomes that contain a SNiP of the sequence (only one varying base)
Functionality
  Load organisms’ genomes from data file
  Returns vector of Genome objects, one for each DNA sequence in the file
  Can return the DNA subsequence starting at some position and is some number of bases long
  Can return name of organism (key)

GenomeMatcher class
Purpose: Given a genome of new organism, identify all genomes with high percentage of matching DNA
Functionality:
  Add new organism’s Genome object to the library
  Find names of all genomes that contain a given sequence or some SNiP, and detail how far from beginning the sequences are
  Given a Genome object, find all genomes that contain at least T% overlap with that obj
    T is determined by user
Must use Trie class to index and search through genomes

GIVEN:
Text files with genomes of common archaea
Simple main driver program to run tests on code
provided.h declares Genome, GenomeMatcher
Both classes are ‘implemented’
All Genome functions are actually implemented in GenomeImpl, which we have to write
All GenomeMatcher functions are to be implemented by us in GenomMatcherImpl

RESTRICTIONS:
No source file other than Genome.cpp may contain name GENOMEIMPL or any helper stuff introduced in Genome.cpp
Can only use the Genome interface provided by SMALLBERG in other class implementations
Trie class template implementations are written directly in header file
Cannot change main.cpp or provided.h

TURNING IN:
Trie.h			      => trie map class template implementation
Genome.cpp		    => Genome class implementation
GenomeMatcher.cpp	=> GenomeMatcher class implementation
report.docx
  Whether any classes have bugs/problems
  Whether or not each method satisfies big O requirements
    If not, what did I do instead and what is new big O
  Pseudocode for two methods
    Trie’s find()
    GenomeMatcher’s findGenomesWithThisDNA()
