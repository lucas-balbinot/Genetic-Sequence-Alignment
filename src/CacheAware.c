#include "Needleman-Wunsch-itmemo.h"

#include <math.h>
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */


/*
 * CACHE DEFINITIONS 
 * they will be set up for the run in the valgrind call
*/
#define SIZE_Z 4096
#define SIZE_L 64

/* EditDistance_CA : It is the main function, performs the calculations
 * in an iterative manner. It dones so by filling up a table that represents
 * the insertions, deletions or substitions, giving the costs of .h
 */
long EditDistance_CA(char* A, size_t lengthA, char* B, size_t lengthB) {

    _init_base_match();

    // make sure A is the longest substring (by changing A and B if necessary)
    if(lengthA<lengthB) {
        char *aux = B;
        B = A;
        A = aux;

        size_t aux_size = lengthA;
        lengthA = lengthB;
        lengthB = aux_size;
    }

    // create the list that will contain the elements
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

    // get the limits (ignoring what is not a base)
    int totalRows = row-1;
    int totalCols = col-1;

    // printf("row: %zu\ncol: %zu\n", totalRows, totalCols);

    // for(int i=0; i<totalRows+1;i++) {
    //     for(int j=0; j<totalCols+1; j++) {
    //         printf("%zu ", edit_dist[i][j]);
    //     }
    //     printf("\n");
    // }

    // divide by the number os bytes that the data type occupies
    // in this case, as char is only 1, this line won't make much effect
    int Z = SIZE_Z / sizeof(char);
    int L = SIZE_L / sizeof(char);

    // initialize the block size to sqrt(Z/2)
    int K = pow(Z/2, 0.5);

    // the 2 outer loops are the blocking parts
    int mid_col, mid_row;

    col = 1;
    for(int i=0; i<lengthA; i+=K) {
        int endA = (i+K < lengthA) ? i+K : lengthA;

        row = 1;
        for(int j=0; j<lengthB; j+=K) {
            int endB = (j+K < lengthB) ? j+K : lengthB;

            mid_col = col;
            // the 2 inner loops are the logic part
            for(int k=i; k<endA; k++) {
                // if A[k] is not a base, no need to enter the loop
                if(!isBase(A[k])) continue;

                mid_row = row;
                for(int l=j; l<endB; l++) {

                    // printf("A[k]:%c | B[l]:%c | i:%d | j:%d | k:%d | l:%d\n",A[k],B[l],i,j,k,l);

                    if(isBase(A[k]) && isBase(B[l])) {
                        // initialization  with cas 1
                        long min = isUnknownBase(A[k]) ?  
                                        SUBSTITUTION_UNKNOWN_COST : 
                                        ( isSameBase(A[k], B[l]) ? 0 : SUBSTITUTION_COST )
                                    + edit_dist[mid_row-1][mid_col-1];
                        {
                            long cas2 = INSERTION_COST + edit_dist[mid_row][mid_col-1] ;      
                            if (cas2 < min) min = cas2 ;
                        }
                        { 
                            long cas3 = INSERTION_COST + edit_dist[mid_row-1][mid_col];      
                            if (cas3 < min) min = cas3 ; 
                        }
                        // printf("  MIN: %d\n", min);
                        edit_dist[mid_row][mid_col] = min;
                        mid_row++;
                    }
                }
                if(isBase(A[k])) mid_col++;
            }
            row = mid_row;
        }
        col = mid_col;
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
