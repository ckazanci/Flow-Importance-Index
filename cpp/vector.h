/* ************************************************************************************************
 * Author: Kelly Black
 *         University of Georgia, Dept. of Mathematics, kjblack@gmail.com
 * Date: October 2018
 *
 * Program to calculate the relative weights associated with each flow in a tropic system
 * with respect to the relative importance. Formula developed with Caner Kazanci, Malcolm
 * Adams, Stuart Whipple, Aladeen Al Basheer, and Bernie Patton.
 *
 * This is a class that is used to keep track of any information organized as either a
 * vector or an array. The SquareMatrix class is used to manipulate the square matrix
 * composed of the columns under consideration. It uses LAPACK to determine if the
 * matrix is of full rank and to approximate the condition number of the matrix.
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
#ifndef VECTOR_H
#define VECTOR_H

#include <ostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <list>

#include "foundfeasible.h"


extern "C" {
    //extern void dswap_(int*,double*,int*,double*,int*);
    //extern void daxpy_(int*,double*,double*,int*,double*,int*);
    //extern void dscal_(int*,double*,double*,int*);
    extern void dgetrf_(int*,int*,double*,int*,int*,int*);
    extern void dgecon_(const char*,int*,double*,int*,double*,double*,double*,int*,int*);
    extern double dlange_(const char*,int*,int*,double*,int*,double*);
}


template <class field>
class Vector
{
public:
    Vector(int number,field initial=0)
    {
        u = new field[number];
        length = number;
        for(int lupe=0;lupe<number;++lupe)
        {
            u[lupe] = initial;
        }
    }

    ~Vector()
    {
      delete [] u;
    }

    field& operator [] (int which)
    {
        return(u[which]);
    }

    field operator [] (int which) const
    {
        return(u[which]);
    }

    friend std::ostream& operator<<(std::ostream& os, Vector& v)
    {
        os << v[0];
        for(int lupe=1;lupe<v.getLength();++lupe)
            os << "," << v[lupe];
        return(os);
    }

    explicit operator field*() {return u;}

    int getLength() { return(length);}

    // Routine to print out a vector to the console
    void printVector()
    {
        std::cout << u[0];
        for(int lupe=1;lupe<getLength();++lupe)
            std::cout << "," << u[lupe];
        std::cout << std::endl;
    }


private:
    field *u;
    int length;
};


template <class field>
class Matrix
{
public:
    // Default constructor
    Matrix()
    {
        u = NULL;
        rows = -1;
        columns = -1;
        index = NULL;
    }

    // Constructor that takes the number rows, columns, and initial value.
    Matrix(int numRows,int numColumns,field initial)
    {
        rows    = numRows;
        columns = numColumns;
        createArray();

        // initialize the values of the entries in the array.
        for(int rowLupe=0;rowLupe<rows;++rowLupe)
        {
            for(int columnLupe=0;columnLupe<columns;++columnLupe)
            {
                u[rowLupe][columnLupe] = initial;
            }
        }
        createIndexPermutation();
    }

    // Constructor that takes the name of a file and reads the
    // matrix from the file.
    Matrix(std::string fileName)
    {
        std::ifstream fp(fileName); // Open a file to read.
        std::string inputLine;
        std::list<std::string> allLines;
        const std::string delimiter = " ";
        std::list<std::string>::iterator oneRow;


        // erase the list of known and unknowable columns.
        unknowableColumns.clearList();
        knownColumns.clearList();

        // Read in the stoichiometry matrix.
        std::getline(fp,inputLine);
        rows = 0;
        while(fp)
        {
            // Add this line to the list of data file lines.
            allLines.push_back(inputLine);
            if(rows<=0)
            {
                // We are still reading the stoichiometry matrix.
                if((inputLine.length()==0)||(inputLine.find(':')!=std::string::npos))
                    rows *= -1; // We have hit the end of the matrix. Stop now.
                else
                    rows -= 1;  // decrement the row count.
            }
            std::getline(fp,inputLine);
        }
        rows = abs(rows);


        // Figure out how many columns there are....
        std::size_t pos = determineNextNondelimiter(0,allLines.front(),delimiter);
        columns = 0;
        while((pos<allLines.front().length())&&(pos != std::string::npos))
        {
            // move past the current column, and get the start of the next column.
            columns += 1;
            pos = allLines.front().find(delimiter,pos+1);

            // At this point we could be at the start of another number,
            // at the end of the string, or at the start of a delimiter.
            // (For example, there may be multiple spaces between numbers.)
            pos = (pos==std::string::npos)?allLines.front().length() :
                          determineNextNondelimiter(pos,allLines.front(),delimiter);
        }
        createArray();

        // Now parse each line to get the entries for each row.
        int currentRow = 0;
        for(oneRow=allLines.begin();(oneRow!=allLines.end());++oneRow)
        {
            if(currentRow<rows)
            {
                // We are still reading the stoichiometry matrix.
                // Parse the line and enter each number.
                pos = determineNextNondelimiter(0,*oneRow,delimiter);
                pos = oneRow->find(delimiter,pos);  // Figure out where the delimiter is.;
                std::size_t currentComma = 0;   // Start reading from the beginning of the string.
                int currentColumn = 0;
                while(pos != std::string::npos)
                {
                    // There is another delimiter. Parse the next number.
                    u[currentRow][currentColumn++] = std::stod(oneRow->substr(currentComma,pos-currentComma));       // Save the number in the array
                    currentComma = pos+1;                                            // update where to start the next search for a comma.
                    currentComma = determineNextNondelimiter(pos,*oneRow,delimiter); // Move past extra delimiters
                    pos = oneRow->find(delimiter,currentComma);                      // Figure out where the comma is.
                }

                // Add the last number in the list to the matrix.
                u[currentRow][currentColumn++] = std::stod(oneRow->substr(currentComma));  // Save the last number in the array
                currentRow += 1;                                                           // move on to the next row.
            }

            else if (oneRow->find("unknowable:") != std::string::npos)
            {
                // This is a row of flows in the stoichiometry matrix that.
                // cannot be measured. Each number coresponds to a column in
                // the stoichiometry matrix whose flow is assumed to be unmeasurable.
                addUnknowableColumns(*oneRow,delimiter);
            }

            else if (oneRow->find("known:") != std::string::npos)
            {
                // This is a row of known flows in the stoichiometry matrix.
                // Each number coresponds to a column in the stoichiometry
                // matrix whose flow is assumed to be already known.
                addKnownColumns(*oneRow,delimiter);
            }
        }

        createIndexPermutation();
        fp.close();
    }

    // Constructor based on a passing another instance of a matrix.
    Matrix(Matrix<field>& other)
    {
        // Set the number of rows and columns.
        rows = other.getNumberRows();
        columns = other.getNumberColumns();
        createArray();

        // Now copy everything in the array over....
        for(int rowLupe=0;rowLupe<rows;++rowLupe)
            for(int colLupe=0;colLupe<columns;++colLupe)
            {
                u[rowLupe][colLupe] = other[rowLupe][colLupe];
            }

        createIndexPermutation();

        // Now copy over the values of the list of unknowables over to this object.
        for(other.beginUnknowableIterations();other.unknowableIterationDone();other.nextUnknowableIteration())
        {
            pushUnknowableValue(other.getCurrentUnknowableValue());
        }

        // Now copy over the values of the list of the known columns over to this object.
        for(other.beginKnownIterations();other.knownIterationDone();other.nextKnownIteration())
        {
            pushKnownValue(other.getCurrentKnownValue());
        }

    }

    // routine to allocate and initialize the space for the array
    void createArray()
    {
        // Allocate the space used by the array.
        u = new field*[rows];
        u[0] = new field[rows*columns];
        for(int rowLupe=0;rowLupe<rows;++rowLupe)
            u[rowLupe] = u[0] + rowLupe*columns;
    }


    // Destructor that deletes the memory allocated in the vector u.
    ~Matrix()
    {
        delete u[0];
        delete [] u;
        rows = -1;
        columns = -1;
    }

    // Method to determine the position of the next character in the string that is not
    // the delimiter
    std::size_t determineNextNondelimiter(std::size_t pos,std::string currentLine,const std::string &delimiter)
    {
        while((pos<currentLine.length())&&(currentLine.substr(pos,1)==delimiter))
        {
            // This position in the string is a delimiter. Check the next character...
            pos += 1;
        }
        return(pos);
    }

    // routine to initialize the index vector that holds the current row permutations.
    void createIndexPermutation()
    {
        // Create the index vector that holds the row permutations.
        index = new int[rows];
        int *ptr = index;
        for(int lupe=0;lupe<rows;++lupe)
            *ptr++ = lupe;  // set the inital value to be the same row number so not permutations.
    }

    void addUnknowableColumns(std::string unknowable,const std::string delimiter)
    {
        std::size_t pos = unknowable.find(":");  // Figure out where the delimiter is.;
        std::size_t currentComma = pos;          // Start reading from the beginning of the string.
        while(pos != std::string::npos)
        {
            // There is another delimiter. Parse the next number.
            if(pos-currentComma>0)
            {
                unknowableColumns.addColumn(std::stod(unknowable.substr(currentComma,pos-currentComma)));
            }
            currentComma = pos+1;                          // update where to start the next search for a comma.
            pos = unknowable.find(delimiter,currentComma); // Figure out where the comma is.
        }

        // Add the last number in the list to the matrix.
        unknowableColumns.addColumn(std::stod(unknowable.substr(currentComma)));
    }

    void addKnownColumns(std::string known,const std::string delimiter)
    {
        std::size_t pos = known.find(":");  // Figure out where the delimiter is.;
        std::size_t currentComma = pos;     // Start reading from the beginning of the string.
        while(pos != std::string::npos)
        {
            // There is another delimiter. Parse the next number.
            if(pos-currentComma>0)
            {
                knownColumns.addColumn(std::stod(known.substr(currentComma,pos-currentComma)));
            }
            currentComma = pos+1;                     // update where to start the next search for a comma.
            pos = known.find(delimiter,currentComma); // Figure out where the comma is.
        }

        // Add the last number in the list to the matrix.
        knownColumns.addColumn(std::stod(known.substr(currentComma)));
    }

    // Method to check to see if any of the entries in the list of unknowables
    // shows up in a given set of indicies.
    bool columnEntryInUnknowable(Vector<int> *indicies,int depth)
    {
        return(columnEntryInList(&unknowableColumns,indicies,depth));
    }

    // Generic method to see if any of the entries in a list shows up
    // in a given set of indices
    bool columnEntryInList(FoundFeasible *columnList,
                           Vector<int> *indicies,int depth)
    {
        // go through each item in the list of indices.
        for(int lupe=0;(lupe<indicies->getLength())&&(lupe<depth);++lupe)
        {
            // if this item exists in the list then return true.
            if(columnList->columnExists((*indicies)[lupe]))
                return(true);
        }
        // looks like the item was not found.
        return(false);
    }

    // Method to check to see if any of the entries in the list of unknowables
    // shows up for one index value
    bool columnEntryInUnknowable(field value)
    {
        // if this item exists in the list then return true.
        return(unknowableColumns.columnExists(value));
    }

    // Method to check to see if any of the entries in the list of knowns
    // shows up for one index value
    bool columnEntryInKnown(field value)
    {
        // if this item exists in the list then return true.
        return(knownColumns.columnExists(value));
    }


    // Method to check to see if *all* of the entries in the list of knowns
    // shows up in a given set of indicies.
    bool allColumnsInKnown(Vector<int> *indicies,int depth)
    {
        return(allColumnsInList(&knownColumns,indicies,depth));
    }

    // Method to check to see if *all* of the entries in the list of knowns
    // shows up in a given set of indicies.
    bool allColumnsInUnknowable(Vector<int> *indicies,int depth)
    {
        return(allColumnsInList(&unknowableColumns,indicies,depth));
    }


    // Generic method to check to see if *all* of the entries in
    // a list shows up in a given set of indicies
    bool allColumnsInList(FoundFeasible *columnList,
                          Vector<int> *indicies,int depth)
    {
        // if there are no entries in the list of known columns
        // default to true.
        if(columnList->length()==0)
            return(true);
        return(columnList->allColumnsExist(indicies,depth));
    }


    // Routine to get a pointer to one row of the matrix using the [] operator
    field*& operator [] (int which)
    {
        if(which >= rows)
        {
            std::cout << "Row number out of bounds" << std::endl;
            exit(1);
        }
        return(u[which]);
    }

    // Routine to get a pointer to one row of the matrix using the [] operator
    field* operator [] (int which) const
    {
        if(which >= rows)
        {
            std::cout << "Row number out of bounds" << std::endl;
            exit(1);
        }
        return(u[which]);
    }

    // Routine used for assignment
    Matrix& operator= (Matrix<field>& other)
    {
        if(this != &other)
        {
            if(rows >= 0)
            {
                // Need to delete the current memory....
                delete u[0];
                delete [] u;
            }

            // Set the number of rows and columns.
            rows = other.getNumberRows();
            columns = other.getNumberColumns();

            // Allocate the space used by the array.
            u = new field*[rows];
            u[0] = new field[rows*columns];
            for(int rowLupe=0;rowLupe<rows;++rowLupe)
                u[rowLupe] = u[0] + rowLupe*columns;

            // Now copy everything in the matrix over....
            for(int rowLupe=0;rowLupe<rows;++rowLupe)
            {
                index[rowLupe] = other.getRowIndex(rowLupe);
                for(int colLupe=0;colLupe<columns;++colLupe)
                {
                    u[rowLupe][colLupe] = other[rowLupe][colLupe];
                }
            }

            // erase the list of known and unknowable columns.
            unknowableColumns.clearList();
            knownColumns.clearList();

            // Now copy over the values of the list of unknowables over to this object.
            for(other.beginUnknowableIterations();other.unknowableIterationDone();other.nextUnknowableIteration())
            {
                pushUnknowableValue(other.getCurrentUnknowableValue());
            }

            // Now copy over the values of the list of the known columns over to this object.
            for(other.beginKnownIterations();other.knownIterationDone();other.nextKnownIteration())
            {
                pushKnownValue(other.getCurrentKnownValue());
            }
        }

        return(*this);
    }

    // Routine to get the data when casting the type
    explicit operator field*() {return u[0];}

    // routine to get the number of rows in the matrix
    int getNumberRows() { return(rows);}

    // routine to get the number of columnss in the matrix
    int getNumberColumns() { return(columns);}

    // Routine to get the value of the index vector at a certain place.
    int getRowIndex(int which)
    {
        return(index[which]);
    }

    /* **************************************************************
     * Routine to swap two given rows in the matrix.
     * ************************************************************** */
    void swapRows(int firstRow,int secondRow)
    {
        // Set up the pointers to point at the first entry in each row.
        field *ptr1 = u[firstRow];
        field *ptr2 = u[secondRow];

        // intermediate values and loop variables required for the loop.
        field tmp;
        int lupe;

        for(lupe=0;lupe<columns;++lupe)
        {
            // go through each entry in the vectors and exchange them.
            tmp = *ptr1;
            *ptr1++ = *ptr2;
            *ptr2++ = tmp;
        }
    }

    /* **************************************************************
     * Routine to add one row of the matrix multipled by a constant
     * to another row in the matrix.
     * ************************************************************** */
    void daxpy(field scaleValue,int changedRow,int sourceRow,int startColumn)
    {
        // Set up the pointers to point at the first entry in each row.
        field *ptr1 = u[changedRow]+startColumn;
        field *ptr2 = u[sourceRow]+startColumn;

        for(int lupe=startColumn;lupe<columns;++lupe)
        {
            // go through each entry and perform the scale/add opteration
            *ptr1++ += scaleValue*(*ptr2++);
        }
    }

    /* **************************************************************
     * Routine to go through a row in the matrix and multiply by a
     * scalar value.
     * ************************************************************** */
    void dscal(field scaleValue,int whichRow,int startColumn)
    {
        // set up the pointer to point at the first entry in the row
        field *ptr = u[whichRow] + startColumn;

        // scale every entry in the vector.
        for(int lupe=startColumn;lupe<columns;++lupe)
            *ptr++ *= scaleValue;
    }

    /* **************************************************************
     * Routine to determine the Reduced Row Echelon Form of the matrix.
     * Performs the RREF in place with the current matrix.
     * ************************************************************** */
    void RREF()
    {
        // define various required variables.
        int outerLupe;
        int innerLupe;
        int currentPivotColumn = 0; // used to indicate the current pivot column

        // Go through each row in the matrix.
        // Figure out a new pivot row and then zero out the entries in the column
        // and and below the current row.
        for(outerLupe=0;(outerLupe<rows) && (currentPivotColumn<columns);++outerLupe)
        {

            // First figure out the current pivot.
            innerLupe = outerLupe; // assume that the pivot is in the first row available.
            while(currentPivotColumn < columns)
            {
                // Start with the current row and work down.
                if((fabs(u[innerLupe][currentPivotColumn])>1e-9))
                {
                    // This entry in this column and row non-zero. Stop here and use this.
                    break;
                }
                else
                {
                    // The entry in this column is essentially zero.
                    innerLupe += 1;    // check the next row.
                    if(innerLupe>=rows)
                    {
                        // We have hit the last row. Everything must be zero. Move over to the next column and start over.
                        innerLupe = outerLupe;
                        currentPivotColumn += 1;
                    }
                }
            }

            if (currentPivotColumn<columns)
            {
                // The current row and pivot are valid.
                // Can zero out the other rows in the current column

                // First see if we need to swap rows.
                if(innerLupe != outerLupe)
                {
                    // The next non-zero row is not the current row.
                    swapRows(innerLupe,outerLupe);
                }

                // zero out everything above the current pivot
                for(int lupe=0;lupe<outerLupe;++lupe)
                {
                    if(fabs(u[lupe][currentPivotColumn]) > 1e-9 )
                    {
                        // Need to zero out the column in this row....
                        daxpy(-u[lupe][currentPivotColumn]/u[outerLupe][currentPivotColumn],
                              lupe,outerLupe,currentPivotColumn);
                    }
                }

                // zero out everything below the current pivot
                for(int lupe=outerLupe+1;lupe<rows;++lupe)
                {
                    if(fabs(u[lupe][currentPivotColumn]) > 1e-9 )
                    {
                        // Need to zero out the column in this column....
                        daxpy(-u[lupe][currentPivotColumn]/u[outerLupe][currentPivotColumn],
                              lupe,outerLupe,currentPivotColumn);
                    }
                }

                // scale the row so that the entry in pivot column is equal to one.
                dscal(1.0/u[outerLupe][currentPivotColumn],outerLupe,currentPivotColumn);

            }

            // About to go on to the next row. Update the current pivot column to use the next column
            currentPivotColumn += 1;
        }


    }


    // Function to print out the whole array
    void printArray()
    {
        int innerLupe;
        int outerLupe;


        std::cout << std::endl << std::endl << getNumberRows() << "-" << getNumberColumns() << std::endl << "     ";

        // Print out the top row that has column numbers
        for(innerLupe=0;innerLupe<getNumberColumns();++innerLupe)
        {
            std::cout <<  "(" << std::setw(2) << innerLupe << ") ";
        }
        std::cout << std::endl;

        for (outerLupe=0;outerLupe<getNumberRows();++outerLupe)
        {
            std::cout << std::endl << "(" << outerLupe << ") " ;
            for(innerLupe=0;innerLupe<getNumberColumns();++innerLupe)
            {
                std::cout <<  std::setw(4) << u[outerLupe][innerLupe]  << " ";
            }
        }
        std::cout << std::endl;

        if(unknowableColumns.length()>0)
        {
            std::cout << "Unknowable: ";
            unknowableColumns.printList();
        }

        if(knownColumns.length()>0)
        {
            std::cout << "Known: ";
            knownColumns.printList();
        }
    }


    // Method to start the process for stepping through the unknowable values.
    void beginUnknowableIterations()
    {
        unknowableColumns.startIteration();
    }

    int getCurrentUnknowableValue()
    {
        return(unknowableColumns.currentValue());
    }

    void nextUnknowableIteration()
    {
        unknowableColumns.next();
    }

    bool unknowableIterationDone()
    {
        return(unknowableColumns.iterationDone());
    }

    void pushUnknowableValue(int value)
    {
        unknowableColumns.addColumn(value);
    }

    // Method to start the process for stepping through the known values.
    void beginKnownIterations()
    {
        knownColumns.startIteration();
    }

    int getCurrentKnownValue()
    {
        return(knownColumns.currentValue());
    }

    void nextKnownIteration()
    {
        knownColumns.next();
    }

    bool knownIterationDone()
    {
        return(knownColumns.iterationDone());
    }

    void pushKnownValue(int value)
    {
        knownColumns.addColumn(value);
    }


protected:
    field **u = NULL;  // the matrix itslef.
    int rows = -1;     // number of rows in the matrix
    int columns = -1;  // number of columns in the matrix
    int* index;        // variable used to hold the current row permutations.

    FoundFeasible unknowableColumns; // list of columns that are considered to be unknowable.
    FoundFeasible knownColumns;      // list of columns that are considered to be known.

};


template <class field>
class SquareMatrix : public Matrix<field>
{
public:
    SquareMatrix() : Matrix<field>()
    {
        work = NULL;
        iwork = NULL;
    }

    SquareMatrix(int numberRows,field initialValue=0) : Matrix<field>(numberRows,numberRows,initialValue)
    {
        createWorkspace();
    }

    SquareMatrix(SquareMatrix<field>& A) : Matrix<field>(A)
    {
        //createIndexPermutation();
        createWorkspace();
    }

    ~SquareMatrix()
    {
        // clean up the mess.
        if(work!=NULL)
            delete [] work;
        if(iwork!=NULL)
            delete [] iwork;
    }

    // method to allocate space for the workspace array
    void createWorkspace()
    {
        work  = new field[5*this->rows];
        iwork = new int[2*this->rows];
    }

    /* *************************************************************************
     * Routine to copy over the columns of a given matrix into the rows
     * of this matrix. The columns to copy are given in a vector of integers.
     * ************************************************************************* */
    void copyColumnsToRows(Matrix<field> *source,Vector<int> *indicies)
    {
        field* ptr;
        for(int rowLupe=0;rowLupe<this->rows;++rowLupe)
        {
            ptr = this->u[rowLupe];
            for(int lupe=0;lupe<this->columns;++lupe)
                *ptr++ = (*source)[rowLupe][(*indicies)[lupe]];
        }
    }

    // http://www.netlib.org/lapack/lug/node80.html
    // http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html#ga5ee879032a8365897c3ba91e3dc8d512
    // http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga188b8d30443d14b1a3f7f8331d87ae60.html#ga188b8d30443d14b1a3f7f8331d87ae60


    /* *************************************************************************
     * Routine to perform an LU decomposition on the current matrix in place.
     * Calls the LAPACK dgetrf_ routine.
     * ************************************************************************* */
    int dgetrf()
    {
        // Set up the variables to pass to the LAPACK routine.
        int numRows = this->rows;
        int numCols = numRows;
        int LDA = numCols;
        int result;

        // Call LAPACK's dgetrf_ routine.
        dgetrf_(&numRows,&numCols,this->u[0],&LDA,this->index,&result);

        // print out the row permutations to see what happened...
        //std::cout << std::endl;
        //for(int lupe=0;lupe<this->rows;++lupe)
        //    std::cout << result << ":" << this->index[lupe] << std::endl;

        // Clean up and return the result.
        return(result);
    }

    /* *************************************************************************
     * Routine to approximate the norm of the current matrix.
     * Calls the LAPACK dlange_ routine.
     * ************************************************************************* */
    double dlange()
    {
        // Set up the parameters to send to the LAPACK method.
        const char *whichNorm = "1";
        int numberRows = this->rows;
        int numberCols = numberRows;
        int LDA = numberRows;

        // Call the lapack routine to calc. the norm of the matrix.
        double norm = dlange_(whichNorm,&numberRows,&numberCols,this->u[0],&LDA,work);

        // Clean up and return the norm.
        return(norm);
    }

    /* *************************************************************************
     * Routine to approximate the reciprocal of the condition number of the
     * current matrix.
     * Calls the LAPACK dgecon_ routine.
     * ************************************************************************* */
    double dgecon(int performLU)
    {
        // see if the lu decomposition has been performed. If necessary
        // go ahead and decompose the matrix.
        if(performLU)
        {
            if(dgetrf()!=0)
                return(0.0);
        }

        // Determine the necessary parameters for the call to the routine
        // to estimate the condition number.
        const char which = '1';
        int numRows = this->rows;
        int LDA = numRows;
        double norm = 1.0;
        double cond = 0.0;
        int result = 0;

        // Call dgecon to get the condition number.
        dgecon_(&which,&numRows,this->u[0],&LDA,&norm,&cond,work,iwork,&result);

        // return the condition number.
        if(result==0)
            return(cond);
        else
            return(0.0);
    }

protected:
    field *work;
    int *iwork;
};

#endif // VECTOR_H
