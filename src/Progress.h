/*    
 *    This file is part of the lidR package by Jean-Romain Roussel et al. and it is 
 *    distributed under the GPL-3 license. If you have questions on its use, please 
 *    or visit:
 *    https://github.com/Jean-Romain/lidR
 */

#ifndef PROGRESS_H
#define PROGRESS_H

#include <Rcpp.h>
#include <time.h>

class Progress
{
  public:
    Progress(unsigned int, std::string = "");
    void check_abort();
    bool check_interrupt();
    void update(unsigned int);
    void increment();
    void exit();
    unsigned int get_iter();

  private:
    unsigned int iter;
    unsigned int iter_max;
    unsigned int percentage;
    unsigned int j;
    std::string prefix;
    clock_t ti;
    clock_t delay;
    bool display;
};

#endif //PROGRESS_H
