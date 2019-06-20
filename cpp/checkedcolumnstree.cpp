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

#include "checkedcolumnstree.h"

// The default constructor. Just assumes there is no tree yet, and the
// width and depth have to be created later.
// At this point this is problematic to use as a constructor.
// There is not a way to correctly set up the root node another way yet,
// so this constructor SHOULD NOT BE USED!
CheckedColumnsTree::CheckedColumnsTree()
{
    setWidth(0);
    setDepth(0);
    checked = nullptr;
}

// Constructor when the max width and depth are known.
// This will create the root node.
CheckedColumnsTree::CheckedColumnsTree(int width,int depth)
{
    setWidth(width);
    setDepth(depth);

    // Now create the root node and initialize each edge to be a null pointer.
    checked = new ColumnTree[width];
    for(int lupe=0;lupe<width;++lupe)
    {
        checked[lupe].value = false;
        checked[lupe].nextColumn = nullptr;
    }
}

// Method used to check to see if a set of columns has been previously checked.
// It also adds the set being checked if it is a new set.
// The order of the columns matters, but the FoundFeasible class automatically
// sorts the list of columns.
bool CheckedColumnsTree::checkColumn(FoundFeasible *indicies)
{
    // Set the number of edges per node and the depth of the tree.
    int width  = getWidth();
    int depth  = getDepth();

    // go through the tree and see if this thing exists. If not add it.
    ColumnTree *current = checked;
    bool exists = true;
    indicies->startIteration();
    for(int depthLupe=0;depthLupe<depth-1;++depthLupe)
    {
        // for each item in indicies see if the level in the tree exists.
        // Note that we know the full depth of the tree in advance which allows
        // to use prev. for loop.
        if(current[indicies->currentValue()].nextColumn==nullptr)
        {
            // This branch of the tree does not exist.
            // Create a new branch and add it to the tree.
            ColumnTree *newBranch = new ColumnTree[width];
            for(int lupe=0;lupe<width;++lupe)
            {
                newBranch[lupe].nextColumn = nullptr;
                newBranch[lupe].value = false;
            }
            current[indicies->currentValue()].nextColumn = newBranch;
            exists = false;
        }
        // Move the pointer down to the next level of the tree.
        current = current[indicies->currentValue()].nextColumn;
        indicies->next();
    }

    // Check the last node in the graph.
    exists = exists && current[indicies->currentValue()].value;
    current[indicies->currentValue()].value = true;
    return(exists);
}

void CheckedColumnsTree::setWidth(int width) {maxWidth=width;}
void CheckedColumnsTree::setDepth(int depth) {maxDepth=depth;}
int CheckedColumnsTree::getWidth(){return(maxWidth);}
int CheckedColumnsTree::getDepth(){return(maxDepth);}
