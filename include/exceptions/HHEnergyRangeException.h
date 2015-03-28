#include <iostream>
#include <stdexcept>
#include <string>

class HHEnergyRangeException: public std::runtime_error
{
public:
 HHEnergyRangeException(const std::string message)
   : std::runtime_error(message) {};
};


