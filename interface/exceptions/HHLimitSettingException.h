#ifndef HHLimitSettingException_
#define HHLimitSettingException_

#include "HHEnergyRangeException.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
  class HHLimitSettingException: public HHEnergyRangeException
{
public:
 HHLimitSettingException(const std::string message)
   : HHEnergyRangeException(message) {};
};
}
#endif
