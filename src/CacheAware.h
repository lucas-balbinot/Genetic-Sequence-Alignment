/**
 * \file CacheAware.h
 * \brief iterative cache aware algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 23/10/2022 
 */

#include "Globals.h" /* have all the cost definitions */

/********************************************************************************
 *  Iterative cache aware algorithm 
 */
/**
 * \fn long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \return :  edit distance between A and B }
 *
 * editDistance_ItMemo is a memoized recursive immplementatioin of Needleman-Wunsch algorithm.
 * It allocates the data structure for memoization table and
 * fills in the memoization table.
 * 
 * If lengthA < lengthB, the sequences A and B are swapped.
 *
 */
long EditDistance_CA(char* A, size_t lengthA, char* B, size_t lengthB);
