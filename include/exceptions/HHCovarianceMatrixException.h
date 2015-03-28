#include <iostream>
#include <stdexcept>
#include <string>

class HHCovarianceMatrixException: public std::runtime_error
{
public:
  HHCovarianceMatrixException(const std::string message)
   : std::runtime_error(message) {};
};


