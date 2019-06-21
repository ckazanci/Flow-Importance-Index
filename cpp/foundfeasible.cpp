/* ************************************************************************************************
 * Author: Kelly Black
 *         University of Georgia, Dept. of Mathematics, kjblack@gmail.com
 * Date: October 2018
 *
 * Program to calculate the relative weights associated with each flow in a tropic system
 * with respect to the relative importance. Formula developed with Caner Kazanci, Malcolm
 * Adams, Stuart Whipple, Aladeen Al Basheer, and Bernie Patton.
 *
 * This is a class that is used to help keep track of one set of vectors. It uses an
 * insertion sort to add new columns sequentially. This way the tree used in the
 * CheckedColumnsTree class will not consider order.
 *
 * The class does not keep track of the column vectors in the stoichiometry matrix.
 * Rather, it keeps track of the column number. So "4" refers to the column vector
 * in the matrix with index 4. (The numbering starts at 0.)
 *
 * Copyright © 2018 Kelly Black
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the “Software”), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * ************************************************************************************************ */

#include <iostream>
#include <iomanip>

#include "foundfeasible.h"
#include "vector.h"

// Base constructor. Nothing to do....
FoundFeasible::FoundFeasible()
{

}

// Method to remove all entries from the list of checked columns.
// This removes everything from the columns list.
void FoundFeasible::clearList()
{
    columns.erase(columns.begin(),columns.end());
    //while(columns.size()>0)
    //    columns.pop_back();
}

// Method to print out all of the items being kept track of in this
// list. Used for debugging.
void FoundFeasible::printList()
{
    std::list<int>::iterator listValues = columns.begin();
    if(listValues!=columns.end())
        std::cout << (*listValues++)+1;
    for(;listValues!=columns.end();++listValues)
    {
        std::cout  << "-" << (*listValues)+1;
    }
    std::cout << std::endl;
}

// Method to add a new column number to the current list of columns being tracked.
// It uses an insertion sort to insure that the column numbers are in ascending
// order.
void FoundFeasible::addColumn(int value)
{
    // Start at the beginning of the list and move through it sequentially.
    std::list<int>::iterator listValues;
    for(listValues=columns.begin();listValues!=columns.end();++listValues)
    {
        if(value < *listValues)
        {
            // Need to insert the value before the current entry in the list.
            columns.insert(listValues,value);
            return;
        }
    }
    // if we get here then the value needs to be appended at the end
    columns.push_back(value);
}

// Method to see if the given value exists in the current list
// of columns being tracked. Returns true if it exists, otherwise
// it is false.
bool FoundFeasible::columnExists(int value)
{
    // go through each item in the list of values.
    std::list<int>::iterator listValues;
    for(listValues=columns.begin();(listValues!=columns.end())&&(*listValues<=value);++listValues)
    {
        if(*listValues==value)
        {
            return(true);      // the item was found!
        }
    }
    // The item was not found.
    return(false);
}

// Method to see if all of the current list of values exist somewhere within a
// given set of indicies.
bool FoundFeasible::allColumnsExist(Vector<int> *indicies, int depth)
{
    //  go through the list of values
    std::list<int>::iterator listValues;
    bool foundIt;
    for(listValues=columns.begin();listValues!=columns.end();++listValues)
    {
        // need to see if the current value is located somewhere in the list.
        // If not then return false.
        foundIt = false;
        for(int lupe=0;(!foundIt) && (lupe<indicies->getLength()) && (lupe<depth);++lupe)
        {
            foundIt = (*listValues==(*indicies)[lupe]);
        }
        if(!foundIt)
            return(false); // this value was missing.
    }

    // everything was found. Huzzah!
    return(true);
}

// Overloaded output operator. Allows the list to be printed
// simply by piping it to an output stream.
std::ostream& operator<<(std::ostream& os, FoundFeasible& v)
{    
    std::list<int>::iterator values;
    for(values=v.columns.begin();values!=v.columns.end();++values)
        os << *values << "-";
    return(os);
}

// Method used to see if a given set of column numbers are located in the list of all
// the sets being tracked. This does a sequential sort through the linked list of
// column sets. IT IS VERY SLOW. It is used in debugging to check other approaches
// to tracking the columns that have been checked. It is slow but is reliable.
bool FoundFeasible::match(Vector<int> *indicies)
{
    // See if the two lists have the same length.
    bool found = static_cast<std::size_t>(indicies->getLength()) == columns.size();

    std::list<int>::iterator listValues;
    for(int columnLupe = 0;found && (columnLupe<indicies->getLength());++columnLupe)
    {
        // for each entry in the vector that was passed, see if we can find a match.
        // Note we assume that each number in the list is distinct.
        bool located = false;
        for(listValues=columns.begin();!located && listValues!=columns.end() && *listValues <= (*indicies)[columnLupe];++listValues)
        {
            // Go through the list of numbers that this data structure is keeping track of.
            // See if any one matches.
            located = located || ((*indicies)[columnLupe]==*listValues);
        }
        found = found && located;
    }
    return(found);
}

// returns the number of column sets being tracked.
int FoundFeasible::length()
{
    return(static_cast<int>(columns.size()));
}

std::list<int>::iterator FoundFeasible::begin(){return(columns.begin());}
std::list<int>::iterator FoundFeasible::end(){return(columns.end());}
void FoundFeasible::startIteration(){currentPos=begin();}
void FoundFeasible::next(){++currentPos;}
bool FoundFeasible::iterationDone(){return(currentPos!=end());}
int FoundFeasible::currentValue(){return(*currentPos);}

