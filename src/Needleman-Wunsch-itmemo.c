#include "Needleman-Wunsch-itmemo.h"

#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;


/* EditDistance_NW_It : It is the main function, performs the calculations
 * in an iterative manner. It dones so by filling up a table that represents
 * the insertions, deletions or substitions, giving the costs of .h
 */
long EditDistance_NW_It(char* A, size_t lengthA, char* B, size_t lengthB) {

    _init_base_match();
    struct NW_MemoContext ctx;
    if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
        {  ctx.X = A ;
        ctx.M = lengthA ;
        ctx.Y = B ;
        ctx.N = lengthB ;
    }
    else
    {  ctx.X = B ;
        ctx.M = lengthB ;
        ctx.Y = A ;
        ctx.N = lengthA ;
    }
    size_t M = ctx.M ;
    size_t N = ctx.N ;
    
    {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
        /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
        * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
        * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
        */ 
        ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
        if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
        for (int i=0; i <= M; ++i) 
        {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
            if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
            for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
        }
    }

    // the evaluation loop
    // starts in 1 to make sense of the table organization
    // as the index 0,0 doesn't represent a char
    for(int i=M; i>=0; i--) {
        for(int j=N; j>=0; j--) {

            long res ;
            char Xi = ctx.X[i] ;
            char Yj = ctx.Y[j] ;
            if (i == ctx.M) /* Reach end of X */
            {  if (j == ctx.N) res = 0;  /* Reach end of Y too */
                else res = (isBase(Yj) ? INSERTION_COST : 0) + ctx.memo[i][j+1];
                //+ EditDistance_NW_RecMemo(c, i, j+1) ;
            }
            else if (j == ctx.N) /* Reach end of Y but not end of X */
            {  res = (isBase(Xi) ? INSERTION_COST : 0) + ctx.memo[i+1][j];
                //+ EditDistance_NW_RecMemo(c, i+1, j) ;
            }
            else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
            {  ManageBaseError( Xi ) ;
                res = ctx.memo[i+1][j];
                //res = EditDistance_Rec_CO(c, i+1, j, K1, K2) ;
            }
            else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
            {  ManageBaseError( Yj ) ;
                res = ctx.memo[i][j+1];
                // EditDistance_Rec_CO(c, i, j+1, K1, K2) ;
            }
            else  
            {  /* Note that stopping conditions (i==M) and (j==N) are already stored in ctx.memo (cf EditDistance_NW_Rec) */ 
                long min = /* initialization  with cas 1*/
                        ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                                : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                        )
                        + ctx.memo[i+1][j+1];
                        //+ EditDistance_NW_RecMemo(c, i+1, j+1) ; 
                { long cas2 = INSERTION_COST + ctx.memo[i+1][j];
                    //+ EditDistance_NW_RecMemo(c, i+1, j) ;      
                if (cas2 < min) min = cas2 ;
                }
                { long cas3 = INSERTION_COST + ctx.memo[i][j+1];
                    //+ EditDistance_NW_RecMemo(c, i, j+1) ;      
                if (cas3 < min) min = cas3 ; 
                }
                res = min ;
            }
            ctx.memo[i][j] = res;
        }
    }

    long res = ctx.memo[0][0];

    { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
    }

   return res;
}
