/* ************************************************************************************************
 * Author: Kelly Black
 *         University of Georgia, Dept. of Mathematics, kjblack@gmail.com
 * Date: October 2018
 *
 * Program to calculate the relative weights associated with each flow in a tropic system
 * with respect to the relative importance. Formula developed with Caner Kazanci, Malcolm
 * Adams, Stuart Whipple, Aladeen Al Basheer, and Bernie Patton.
 *
 * This is a class that is used to help keep track of which sets of columns have
 * already been checked. It organizes the set of checked vectors in a tree. When
 * a set of vectors has been checked it also adds the vector to the current
 * tree.
 *
 * The order of the vectors matters, and this will differentiate the columns checked
 * by the order they are given. If you do not want it to consider order, then you should
 * give it the set of vectors in a sorted order.
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
#ifndef CHECKEDCOLUMNSTREE_H
#define CHECKEDCOLUMNSTREE_H

#include <vector>

#include "vector.h"
#include "foundfeasible.h"

class CheckedColumnsTree
{
public:
    CheckedColumnsTree();
    CheckedColumnsTree(int width,int depth);

    bool checkColumn(FoundFeasible *indicies);

    void setWidth(int width);
    void setDepth(int depth);
    int getWidth();
    int getDepth();

    // The information is kept as a tree. The primary data structure
    // used is the following struct. Each node of the tree keeps
    // an array of items of the ColumnTree type. Each node then
    // includes an array of pointers through the nextColumn pointer
    // in the struct.
    struct ColumnTree
    {
        bool value;
        ColumnTree *nextColumn;
    };

private:
    int maxWidth;          // the width of the tree. (Number of edges for each node.)
    int maxDepth;          // the depth of the tree.
    ColumnTree *checked;   // A pointer to the tree's root.
};

#endif // CHECKEDCOLUMNSTREE_H
