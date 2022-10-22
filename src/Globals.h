#include <stdlib.h> /* for size_t */

/*
 * Costs for operations on canonical bases
 * Three  operations: insertion and sustitution of one base by an another 
 * Note= substitution of an unknown base N by another one (known or unknown) as the same cost than substitution between 2 different known bases
 */
/** \def SUBSTITUTION_COST
 *  \brief Cost of substitution of one canonical base by another
 */
#define SUBSTITUTION_COST	1

/** \def SUBSTITUTION_UNKNOWN_COST
 *  \brief Cost of substitution of an unknown base (N) by another one (canonical or unknown)
 */
#define SUBSTITUTION_UNKNOWN_COST	1  /* Cost for sustitition of an Unknown bas N by another on -known or unkown- */ 

/** \def INSERTION_COST
 *  \brief Cost of insertion of a canonical base 
 */
#define INSERTION_COST		2
