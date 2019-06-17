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
#ifndef FOUNDFEASIBLE_H
#define FOUNDFEASIBLE_H

#include <list>
//#include <vector>
//#include "vector.h"

template<class field> class Vector;

class FoundFeasible
{
public:
    FoundFeasible();

    friend std::ostream& operator<<(std::ostream& os, FoundFeasible& v);
    void clearList();
    void printList();
    void addColumn(int value);
    bool columnExists(int value);
    bool allColumnsExist(Vector<int> *indicies,int depth);
    bool match(Vector<int> *indicies);
    int length();

    // Define the iterators used to step through the columns.
    std::list<int>::iterator begin();
    std::list<int>::iterator end();
    void startIteration();
    void next();
    bool iterationDone();
    int currentValue();

private:
    std::list<int> columns;
    std::list<int>::iterator currentPos;
};

#endif // FOUNDFEASIBLE_H
