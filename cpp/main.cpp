/* ************************************************************************************************
 * Author: Kelly Black
 *         University of Georgia, Dept. of Mathematics, kjblack@gmail.com
 * Date: October 2018
 *
 * Program to calculate the relative weights associated with each flow in a tropic system
 * with respect to the relative importance. Formula developed with Caner Kazanci, Malcolm
 * Adams, Stuart Whipple, Aladeen Al Basheer, and Bernie Patton.
 *
 * This is the main program. It reads in a file from the command line that has
 * the stoichiometry matrix. It then determines the RREF of the matrix and goes
 * through all combinations of the column vectors that might form a set of full
 * rank. It then keeps track of which vectors were in the set and maintains the
 * relevant statistics.
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
#include <string>
#include <list>

#include "vector.h"
#include "foundfeasible.h"
#include "checkedcolumnstree.h"

//using namespace std;
//#define DEBUG

// Routine to see if a given column is already listed in the
// vector of indicies to check.
bool columnExists(Vector<int> *indicies,int currentRow,int value)
{
    // Go through all of the previous entries in the vector.
    for(int lupe=0;lupe<currentRow;++lupe)
        if((*indicies)[lupe]==value)  // If the value exists return true. All done!
            return(true);

    // If we get here the value was not found.
    return(false);
}

// Routine to see if a given collection of columns has already been considered.
// The idea is that if the coresponding entries in the RREF are non-zero, and
// the column numbers are in descending order, then the set has already been
// considered.
bool columnsConsidered(Matrix<double> *rref,
                       Vector<int> *indicies,
                       int currentRow,
                       int currentColumn)
{

    bool considered = false;

    //std::cout << "checking: " << (*indicies)[currentRow] << "/" << currentRow << "/" << currentColumn << " " << std::endl;
    //printVector(*indicies);
    for(int prevColumn=0;!considered && (prevColumn<currentRow);++prevColumn)
    {
        if((*indicies)[prevColumn]>(*indicies)[currentRow])
        {
            // This column may have been considered in a previous iteration.
            for(int rowCheck=0;!considered && (rowCheck<indicies->getLength());++rowCheck)
            {
#ifdef DEBUG
                std::cout << std::setw(5) << (*rref)[rowCheck][(*indicies)[rowCheck]]
                          << std::setw(5) << (*rref)[currentRow][(*indicies)[rowCheck]]
                          << std::setw(5) << (*rref)[rowCheck][currentColumn]
                          << std::setw(5) << (*rref)[currentRow][currentColumn]
                               << std::endl;
#endif
                if(rowCheck!=prevColumn)
                     considered = considered ||
                             (fabs((*rref)[prevColumn][(*indicies)[prevColumn]] *
                                (*rref)[prevColumn][currentColumn] *
                                (*rref)[rowCheck][(*indicies)[prevColumn]] *
                                (*rref)[rowCheck][currentColumn]) > 1e-4 );
            }
        }
    }

    // if we get here then the set of columns have not been considered previously.
    return(considered);
}

// http://www.netlib.org/lapack/lug/node80.html
// http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html#ga5ee879032a8365897c3ba91e3dc8d512
// http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga188b8d30443d14b1a3f7f8331d87ae60.html#ga188b8d30443d14b1a3f7f8331d87ae60


/* *******************************************************************************
 * Routine to check to see if a given set of columns has been previously checked.
 * ******************************************************************************* */
bool columnsPreviouslyChecked(Vector<int> *indicies,
                              std::list<FoundFeasible*> *checkedSets)
{
    bool alreadyChecked = false;
    std::list<FoundFeasible*>::iterator prevChecked;
    for(prevChecked=checkedSets->begin();!alreadyChecked && (prevChecked!=checkedSets->end());++prevChecked)
    {
        // go through each list in the set of previously checked columns. If there is a match then stop the show.
        alreadyChecked = alreadyChecked || (*prevChecked)->match(indicies);
    }
    return(alreadyChecked);
}

/* *******************************************************************************
 * Routine to check a given set of columns and get the condition number for
 * the columns from the original matrix.
 * ******************************************************************************* */
bool testFullColumnSet(Matrix<double> *originalStoichiometry,
                       SquareMatrix<double> *testBasis,
                       int *numberFeasible,
                       Vector<int> *feasibleByColumn,
                       std::list<double> *conditionNumbers,
                       Vector<double> *sumConditionNumbers,
                       Vector<double> *sumInvConditionNumbers,
                       Vector<int> *indicies,
                       int *numberRepeats,
                       CheckedColumnsTree   *previouslyChecked)
{
    // form a matrix with the appropriate columns of the original stoichiometry matrix
    // create the list of columns and keep them in a useable form.
    // First record the current combination.
    FoundFeasible newColumns;
    newColumns.clearList();
    for(int foundLupe=0;foundLupe<indicies->getLength();++foundLupe)
    {
        newColumns.addColumn((*indicies)[foundLupe]);
    }

    // See if this combination has been previously checked.
    if(previouslyChecked->checkColumn(&newColumns))
    {
        // We have seen this before. Time to go....
        (*numberRepeats)++;
        return(false);
    }

    // see if the vectors have full rank.
    testBasis->copyColumnsToRows(originalStoichiometry,indicies); // copy the columns in question to a separate matrix.
    SquareMatrix<double> originalCopy(*testBasis);                // make a copy for determining the norm of the matrix later.
    if(testBasis->dgetrf()==0)
    {
        // The resulting system is of full rank.
        // Figure out the necessary details.
        // First get the matrix norm.
        double matrixNorm  = originalCopy.dlange();
        double inverseNorm = testBasis->dgecon(false);
        double conditionNumber  = matrixNorm*inverseNorm;
        conditionNumbers->push_back(conditionNumber);
        //std::cout << *numberFeasible << " " << matrixNorm << "/" << inverseNorm << std::endl;

        //else
        {
            *numberFeasible += 1;    // increment the number of feasible sets.

            // Now increment the values of the feasible columns for each column not
            // given in the current list of indices. Takes advantage that the
            // newColumns variable is sorted as the items were added above.
            int lupe=0;
            for(newColumns.startIteration();newColumns.iterationDone();newColumns.next())
            {
                // go through each column in the current list.
                while((lupe<feasibleByColumn->getLength())&&(lupe<newColumns.currentValue()))
                {
                    // The column number to test is less than the current column number in the list
                    // of columns used to invert the matrix above.
                    (*feasibleByColumn)[lupe] += 1;
                    (*sumInvConditionNumbers)[lupe] += conditionNumber;
                    (*sumConditionNumbers)[lupe] += 1.0/conditionNumber;
                    lupe += 1;
                }
                // We are either at the end of the list or found a column number equal to
                // one in the list to test. Move on to the next column.
                lupe += 1;
            }

            // We have gone through everything, but there might be a few extra columns
            // at the end whose column numbers were larger than what was in the test list.
            // Increment everything that follows.
            while(lupe<feasibleByColumn->getLength())
            {
                (*feasibleByColumn)[lupe] += 1;
                (*sumInvConditionNumbers)[lupe] += conditionNumber;
                (*sumConditionNumbers)[lupe] += 1.0/conditionNumber;
                lupe += 1;
            }

        }
    }

    return(true);
}

/* *******************************************************************************
 * Routine to go through the RREF of the matrix and get all combinations of the
 * columns that have non-zero entries in the RREF of the matrix.
 *
 * This is a recursive routine. If starts with the top row and goes through each
 * column with a non-zero entry. For each of those entries it then calls the routine
 * to go through the next row and check each column for the next row. (Repeat!)
 * ******************************************************************************* */
void checkColumns(Matrix<double> *rref,
                  Matrix<double> *originalStoichiometry,
                  SquareMatrix<double> *testBasis,
                  Vector<int> *indicies,
                  int *numberFeasible,
                  int currentRow,
                  int *numberRepeats,
                  Vector<int> *feasibleByColumn,
                  std::list<double> *conditionNumbers,
                  Vector<double> *sumConditionNumbers,
                  Vector<double> *sumInvConditionNumbers,
                  CheckedColumnsTree   *previouslyChecked
                  )
{
    // Go through each column for the current row. Check to see which entries are
    // non-zero....
    for(int lupe=0;lupe<rref->getNumberColumns();++lupe)
    {
        (*indicies)[currentRow] = lupe;
        if((fabs((*rref)[currentRow][lupe])>1.0e-9)
                && (!columnExists(indicies,currentRow,lupe))
                && (!originalStoichiometry->columnEntryInKnown(lupe))
                //&& (!columnsConsidered(rref,indicies,currentRow,lupe))
                )
        {
            // This entry in the RREF matrix is non-zer0.
            // The column has also not been checked previously.
            // It is a potential column to check.

            // Add this column to the list of columns currently under consideration.
            if( (currentRow+1) >= rref->getNumberRows())
            {
                if(originalStoichiometry->allColumnsInUnknowable(indicies,currentRow+1))
                {
                    // We now have a full list of columns to check.
                    testFullColumnSet(originalStoichiometry,testBasis,
                                      numberFeasible,feasibleByColumn,conditionNumbers,
                                      sumConditionNumbers,sumInvConditionNumbers,
                                      indicies,numberRepeats,previouslyChecked);
                }
            }
            else
            {
                // We need at least one more column to check.
                // search using the next row in the RREF matrix.
                checkColumns(rref,originalStoichiometry,testBasis,indicies,
                             numberFeasible,currentRow+1,numberRepeats,
                             feasibleByColumn,conditionNumbers,sumConditionNumbers,sumInvConditionNumbers,previouslyChecked);
            }
        }
    }
}

/* *******************************************************************************
 * Routine to calculate the value of n choose k. Uses the identity that
 * n choose k is equal to (n/k)*( (n-1) choose (k-1) ).
 * ******************************************************************************* */
long combinations(int n,int k)
{
    // make sure everything is okay.
    if(n<k)
        return(0);
    else if (2*k < n)
        // use the identity that nCk = nCn-k.
        // In this case there are fewer computations for nCn-k.
        return(combinations(n,n-k));

    long combinations = 1;
    long lupe;
    long denominator = 1;
    for(lupe=static_cast<long>(n);lupe>static_cast<long>(n-k);--lupe)
        combinations = combinations*lupe/(denominator++);
    return(combinations);
}



int main(int argc,char **argv)
{
    if(argc < 2)
    {
        std::cout << "Error - Command line should include name of file that has the stoichiometry matrix: \"" << argv[0] << " stoich.txt.\"" << std::endl << std::endl;
        exit(1);
    }

    // Variables used to track the stoich. matrix and the columns being tested
    Matrix<double>       stoichiometry(argv[1]);                          // The stoichiometry matrix read from a file.
    Matrix<double>       originalStoich(stoichiometry);                   // A copy of the stoich. matrix. Used to construct test matrices from its columns.
    SquareMatrix<double> testBasis(stoichiometry.getNumberRows(),0.0);    // A square matrix made up of columns of the stoich. matrix. For testing.
    Vector<int>          indicies(stoichiometry.getNumberRows(),-1);      // List of columns to test to see if they form a set of full rank.

    // Variables used to track the statistics associated with the columns from the stoich. matrix that have been tested.
    Vector<int>          feasibleByColumn(stoichiometry.getNumberColumns(),0);           // Number of times a given column appears as a feasible set.
    Vector<double>       sumConditionNumbers(stoichiometry.getNumberColumns(),0.0);      // Sum of the cond. #'s for each matrix a column appears in a feasible set.
    Vector<double>       sumInvConditionNumbers(stoichiometry.getNumberColumns(),0.0);   // Sum of 1/cond. #'s for each matrix a column appears in a feasible set.
    CheckedColumnsTree   *previouslyChecked =                                            // tree structure uses to track all combinations of columns that have been tested.
            new CheckedColumnsTree(stoichiometry.getNumberColumns(),stoichiometry.getNumberRows());
    std::list<double>    conditionNumbers;                                               // List of all condition numbers for feasible sets.

    // Test variables used for debugging.
    int numberFeasible = 0;   // Total number of feasible sets tested.
    int numberRepeats = 0;    // Number of times vectors were tested that have already been in a feasible test.


    stoichiometry.printArray();
    stoichiometry.RREF();
    //stoichiometry.printArray();
    checkColumns(&stoichiometry,&originalStoich,&testBasis,&indicies,
                 &numberFeasible,0,&numberRepeats,
                 &feasibleByColumn,&conditionNumbers,
                 &sumConditionNumbers,&sumInvConditionNumbers,
                 previouslyChecked);

    long numberPossible = combinations(stoichiometry.getNumberColumns()-1,stoichiometry.getNumberRows());
    std::cout << "Number Feasible: " << numberFeasible << std::endl
              << "Normalization: "   << numberPossible << std::endl
              << "Feasible by column: " << std::endl
              << "Node Feasible     Sum Cond.   Sum Inv Cond         Impact" << std::endl;
    for(int lupe=0;lupe<feasibleByColumn.getLength();++lupe)
    {
        std::cout << std::fixed
                  << std::setw(4) << lupe << "    "
                  << std::setw(5) << feasibleByColumn[lupe] << "   "
                  << std::setw(11) << std::setprecision(5) << sumConditionNumbers[lupe] << "    "
                  << std::setw(11) << std::setprecision(5) << sumInvConditionNumbers[lupe] << "    ";
        if(feasibleByColumn[lupe]>0)
        {
            // This flow appears in at least one valid representation.
            std::cout << std::setw(11) << std::setprecision(5)
                      << sumInvConditionNumbers[lupe]*sumConditionNumbers[lupe]/(static_cast<double>(feasibleByColumn[lupe])*static_cast<double>(numberFeasible));
        }
        else
        {
            // This flow does not appear in any valid representations.
            std::cout << "         NA";
        }
        std::cout<< std::endl;
    }

#ifdef DEBUG
    std::cout << "Number of repeats: " << numberRepeats << std::endl;
#endif

    return 0;
}

