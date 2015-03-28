#include <iostream>
#include <stdexcept>
#include <string>

class HHEnergyConstraintException: public std::runtime_error
{
public:
  HHEnergyConstraintException(const std::string message)
   : std::runtime_error(message) {};
};


