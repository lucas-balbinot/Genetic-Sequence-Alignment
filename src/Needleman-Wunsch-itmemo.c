/*
 * POSSIBLE IMPROVEMENT : the inner loop being the smaller array
 */
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

    _init_base_match();

    // creates the tables of length
    // the +1 to take deletions and insertions into account
    // simulates long edit_dist[lengthA+1][lengthB+1] = {0};
    long** edit_dist = (long**)malloc( (lengthA+1) * sizeof(long*));
    for(int i=0; i<lengthB; i++) {
        edit_dist[i] = (long*)malloc( (lengthB+1) * sizeof(long));
        edit_dist[i] = 0;
    }

    // rows and cols to ignore the chars that are not bases
    int row = 1, col = 1;
    // just insertions and deletions for the first row and column
    for(int i=1; i<lengthA+1; i++) {

        // for the row
        if(isBase(A[i])) {
            // if it is a base, calculate the cost of insertion
            edit_dist[0][col] = INSERTION_COST * i;
            col++;
            // otherwise, don't add anything to the value
        }

    }
    for(int i=1; i<lengthB+1; i++) {
        // for the column
        if(isBase(B[i])) {
            edit_dist[row][0] = INSERTION_COST * i;
            row++;
        }
    }

    // the evaluation loop
    // starts in 1 to make sense of the table organization
    // as the index 0,0 doesn't represent a char
    long delta;
    row = 1;
    for(int i=1; i<lengthA+1; i++) {
        col = 1;
        for(int j=1; j<lengthB+1; j++) {
            // test if it is indeed a base
            if(isBase(A[i]) && isBase(B[j])) {
                // initialization  with cas 1
                long min = isUnknownBase(A[i]) ?  
                                SUBSTITUTION_UNKNOWN_COST : 
                                ( isSameBase(A[i], B[i]) ? 0 : SUBSTITUTION_COST )
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

                // update column count
                col++;
            }
            // update row count
            row++;
        }
    }

    return edit_dist[row][col];

}
