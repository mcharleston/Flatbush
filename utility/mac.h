/*
 * mac.h
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#ifndef MAC_H_
#define MAC_H_

//#define DEBUGGING

#ifdef DEBUGGING
#define DEBUG(A) { if (_debugging) { A; } }
#endif
#ifndef DEBUGGING
#define DEBUG(A) {}
#endif


#endif /* MAC_H_ */
