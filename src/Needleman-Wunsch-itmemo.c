#include "Needleman-Wunsch-itmemo.h"

#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/* EditDistance_NW_It : It is the main function, performs the calculations
 * in an interative manner. It dones so by filling up a table that represents
 * the insertions, deletions or substitions, giving the costs of .h
 */
long EditDistance_NW_It(char* A, size_t lengthA, char* B, size_t lengthB) {

    if(lengthA<lengthB) {
        char *aux = B;
        B = A;
        A = aux;

        size_t aux_size = lengthA;
        lengthA = lengthB;
        lengthB = aux_size;
    }

    _init_base_match();

    // creates the tables of length
    // the +1 to take deletions and insertions into account
    // simulates long edit_dist[lengthA+1][lengthB+1] = {0};

    // lengthA + 1 = number of rows
    long** edit_dist = (long**)malloc( (lengthB+1) * sizeof(long*));
    for(int i=0; i<lengthA+1; i++) {
        // lengthB + 1 = number os columns
        edit_dist[i] = (long*)malloc( (lengthA+1) * sizeof(long));
        for(int j=0; j<lengthA+1; j++)
            edit_dist[i][j] = 0;
    }

    //rows and cols to ignore the chars that are not bases
    int row = 1, col = 1;
    // just insertions and deletions for the first row and column
    for(int i=1; i<lengthA+1; i++) {
        // for the row
        if(isBase(A[i-1])) {
            // if it is a base, calculate the cost of insertion
            edit_dist[0][col] = INSERTION_COST * col;
            // printf("col[%d] -> %zu\n", col, edit_dist[0][col]);
            col++;
            // otherwise, don't add anything to the value
        }
    }
    for(int i=1; i<lengthB+1; i++) {
        // for the column
        if(isBase(B[i-1])) {
            //printf("%c", B[i-1]);
            edit_dist[row][0] = INSERTION_COST * row;
            // printf("row[%d] -> %zu\n", row, edit_dist[row][0]);
            row++;
        }
    }

    int totalCols = col - 1;
    int totalRows = row - 1;

    // printf("row: %zu\ncol: %zu\n", totalRows, totalCols);

    // for(int i=0; i<totalRows+1;i++) {
    //     for(int j=0; j<totalCols+1; j++) {
    //         printf("%zu ", edit_dist[i][j]);
    //     }
    //     printf("\n");
    // }


    // the evaluation loop
    // starts in 1 to make sense of the table organization
    // as the index 0,0 doesn't represent a char
    long delta;
    col = 1;
    for(int i=0; i<lengthA; i++) {
        row = 1;
        for(int j=0; j<lengthB; j++) {
            // test if it is indeed a base
            if(isBase(A[i]) && isBase(B[j])) {
                // initialization  with cas 1
                long min = isUnknownBase(A[i]) ?  
                                SUBSTITUTION_UNKNOWN_COST : 
                                ( isSameBase(A[i], B[j]) ? 0 : SUBSTITUTION_COST )
                            + edit_dist[row-1][col-1];
                { 
                    long cas2 = INSERTION_COST + edit_dist[row][col-1] ;      
                    if (cas2 < min) min = cas2 ;
                }
                { 
                    long cas3 = INSERTION_COST + edit_dist[row-1][col];      
                    if (cas3 < min) min = cas3 ; 
                }
                // the value is updated with the min
                edit_dist[row][col] = min;
                
                // update row count
                row++;
            }
        }
        if(isBase(A[i])) {
            // update column count
            col++;
        }
    }

    // printf("row: %zu\ncol: %zu\n", totalRows, totalCols);

    // for(int i=0; i<totalRows+1;i++) {
    //     for(int j=0; j<totalCols+1; j++) {
    //         printf("%zu ", edit_dist[i][j]);
    //     }
    //     printf("\n");
    // }

    return edit_dist[totalRows][totalCols];
}
