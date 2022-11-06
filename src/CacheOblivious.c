/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


// #include "Needleman-Wunsch-recmemo.h"
#include "CacheOblivious.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <math.h>
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/

#define Z 4096
#define S 30

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

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_Rec_CO(struct NW_MemoContext *c, int begin_1, int begin_2, int end_1, int end_2) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
    int n_1 = end_1 - begin_1;
    int n_2 = end_2 - begin_2;

    if((n_1<S) && (n_2<S)) {

        // printf("begin1:%d begin2:%d\n", begin_1, begin_2);

        for(int i=end_1; i>=begin_1; i--) {    
            for(int j=end_2; j>=begin_2; j--) {

                // printf("i:%02zu | j:%02zu | I:%02zu | J:%02zu | K1:%02zu | K2:%02zu\n", i,j,begin_1,begin_2,end_1,end_2);

                long res ;
                char Xi = c->X[i] ;
                char Yj = c->Y[j] ;
                if (i == c->M) /* Reach end of X */
                {  if (j == c->N) res = 0;  /* Reach end of Y too */
                    else res = (isBase(Yj) ? INSERTION_COST : 0) + c->memo[i][j+1];
                    //+ EditDistance_NW_RecMemo(c, i, j+1) ;
                }
                else if (j == c->N) /* Reach end of Y but not end of X */
                {  res = (isBase(Xi) ? INSERTION_COST : 0) + c->memo[i+1][j];
                    //+ EditDistance_NW_RecMemo(c, i+1, j) ;
                }
                else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
                {  ManageBaseError( Xi ) ;
                    res = c->memo[i+1][j];
                    //res = EditDistance_Rec_CO(c, i+1, j, K1, K2) ;
                }
                else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
                {  ManageBaseError( Yj ) ;
                    res = c->memo[i][j+1];
                    // EditDistance_Rec_CO(c, i, j+1, K1, K2) ;
                }
                else  
                {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
                    long min = /* initialization  with cas 1*/
                            ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                                    : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                            )
                            + c->memo[i+1][j+1];
                            //+ EditDistance_NW_RecMemo(c, i+1, j+1) ; 
                    { long cas2 = INSERTION_COST + c->memo[i+1][j];
                        //+ EditDistance_NW_RecMemo(c, i+1, j) ;      
                    if (cas2 < min) min = cas2 ;
                    }
                    { long cas3 = INSERTION_COST + c->memo[i][j+1];
                        //+ EditDistance_NW_RecMemo(c, i, j+1) ;      
                    if (cas3 < min) min = cas3 ; 
                    }
                    res = min ;
                }
                c->memo[i][j] = res;
            }
        }
    }
    else {
        if(n_1>n_2){
            // printf("N1>N2 -> begin:%zu, end:%zu\n",begin_1,end_1);
            EditDistance_Rec_CO(c, (begin_1+end_1)/2, begin_2, end_1, end_2);
            EditDistance_Rec_CO(c, begin_1, begin_2, (begin_1+end_1)/2, end_2);
        }
        else {
            // printf("N2>N1 -> begin:%zu, end:%zu\n",begin_2,end_2);
            EditDistance_Rec_CO(c, begin_1, (begin_2+end_2)/2, end_1, end_2);
            EditDistance_Rec_CO(c, begin_1, begin_2,end_1, (begin_2+end_2)/2);
        }
    }

    return c->memo[0][0];
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_CO(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
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

   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_Rec_CO( &ctx, 0, 0, M, N ) ;
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

