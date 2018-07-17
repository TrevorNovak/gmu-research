#ifndef HYPERBOX_H_
#define HYPERBOX_H_

#include "masterheader.h"
#include <map>

void hyperbox(bool stocsim, const std::vector<Solution*>& visited, long starIndex, long oldStarIndex, std::vector<double*>& Ak, \
			std::vector<double>& bk, const std::vector<double*>& A, const std::vector<double>& b, \
			std::vector<long>& activeSolutions, const std::vector<CornerMap*>& corners, bool boxalone = false);



#endif /*HYPERBOX_H_*/